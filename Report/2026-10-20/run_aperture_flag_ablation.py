"""Rebuild COSMOS catalogs after dropping only the aperture flag == 0 cut.

Raw SourceXtractor++ catalogs and production results are read-only inputs.
Run field jobs sequentially, then combine, plot and check input fingerprints:
  python run_aperture_flag_ablation.py field 1  # repeat for 2, 3, 4
  python run_aperture_flag_ablation.py combined
  python run_aperture_flag_ablation.py verify

Notebook algorithms are executed directly, rather than reimplemented here.
Only output paths, field number and SDSS cache reuse are execution overrides.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT = ROOT/'output/experiments/aperture_flag_ablation'
BASELINE = EXPERIMENT/'baseline'


def notebook(name):
    return json.loads((ROOT/'Notebooks'/name).read_text())


def cell_source(nb, index):
    return ''.join(nb['cells'][index]['source'])


def execute(nb, indices, namespace):
    for index in indices:
        print(f'CELL {index} START {time.strftime("%H:%M:%S")}', flush=True)
        exec(compile(cell_source(nb, index), f'{namespace["__name__"]}:cell{index}', 'exec'), namespace)
        plt.close('all')
        print(f'CELL {index} DONE', flush=True)


def input_fingerprints():
    """Hash all reused raw catalogs and footprint masks to detect any rewrite."""
    paths = sorted((ROOT/'output').glob('COSMOS*/singlephot/cat_cosmos*_m*.fits'))
    paths += sorted((ROOT/'COSMOS_data').glob('COSMOS*/COSMOS*_validmask.fits'))
    result = {}
    for path in paths:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024*1024), b''):
                digest.update(chunk)
        result[str(path.relative_to(ROOT))] = dict(sha256=digest.hexdigest(), size=path.stat().st_size,
                                                   mtime_ns=path.stat().st_mtime_ns)
    return result


def forbid_photometry(*args, **kwargs):
    """Fail closed if any notebook unexpectedly attempts a shell process."""
    raise RuntimeError('Subprocess execution is disabled: reuse existing raw photometry.')


def run_field(field):
    nb = notebook('2_singphotandmatching.ipynb')
    ns = {'__name__': f'cosmos{field}_flag_ablation'}
    # The setup still reads the exact production images, weights and masks.
    setup = cell_source(nb, 2).replace('COSMOS = 1\n', f'COSMOS = {field}\n', 1)
    exec(compile(setup, 'singlephot:setup', 'exec'), ns)
    product = EXPERIMENT/f'COSMOS{field}/singlephot'
    diagnostic, analysis = product/'diagnostics', product/'analysis'
    diagnostic.mkdir(parents=True, exist_ok=True)
    analysis.mkdir(parents=True, exist_ok=True)
    ns.update(PRODUCT_DIR=product, MATCH_DIR=diagnostic, ANALYSIS_DIR=analysis, MATCH_FIG_DIR=diagnostic)
    execute(nb, [4], ns)
    ns['OVERWRITE_PHOTOMETRY'] = False
    ns['OVERWRITE_SDSS'] = False
    # Cell 6 is deliberately not executed at all; even a missing raw file
    # must fail rather than invoke SourceXtractor++.
    subprocess.run = forbid_photometry
    for band in ns['ACTIVE_BANDS']:
        assert ns['get_catalog_path'](field, band).is_file()
    cache_name = f'cosmos{field}_master_sdss_matching_primary_dr17.fits'
    cache = ROOT/f'output/COSMOS{field}/singlephot/diagnostics'/cache_name
    if cache.exists():
        shutil.copy2(cache, diagnostic/cache_name)
    execute(nb, [8], ns)
    rows = []
    for band in ns['ACTIVE_BANDS']:
        raw, chosen = ns['catalogs'][band], ns['matching_catalogs'][band]
        # The old selection is the new selection intersected with flag == 0;
        # all other conditions are identical, including the valid-mask lookup.
        flags = np.asarray(chosen['aper3_flags'])
        rows.append(dict(COSMOS=field, band=band, n_raw=len(raw),
                         n_selected_baseline=int(np.sum(flags == 0)),
                         n_selected_no_flag=len(chosen), n_added=int(np.sum(flags != 0))))
    Table(rows=rows).write(product/'selection_comparison.csv', format='ascii.csv', overwrite=True)
    execute(nb, [10,12,14,16], ns)
    sdss_code = cell_source(nb, 18)
    definitions, separator, attachment = sdss_code.partition('\nsdss_primary = ')
    assert separator
    exec(compile(definitions, 'singlephot:sdss_definitions', 'exec'), ns)
    original_loader = ns['load_sdss_primary_catalog']

    def reuse_complete_sdss_cache(reference_catalog):
        """Reuse cached candidates only when the complete match halo is covered.

        The notebook downloads a 5-arcsec margin but matches at 1.5 arcsec.
        Tiny new centroid extrema can exceed the original *5-arcsec* bounding
        box while its cached sky region still contains every possible match.
        Prove containment at the unchanged matching radius before reuse;
        never relabel the original cache bounds or change the match radius.
        """
        path=ns['SDSS_PRIMARY_PATH']
        if path.exists():
            cached=Table.read(path)
            ra=np.asarray(reference_catalog['master_ra'],float)
            dec=np.asarray(reference_catalog['master_dec'],float)
            margin=ns['SDSS_MATCH_RADIUS_ARCSEC']/3600
            # Enlarge RA by cos(dec) for a conservative spherical halo.
            ra_margin=margin/np.min(np.cos(np.deg2rad(dec)))
            complete=(float(cached.meta.get('RA_MIN',np.inf)) <= ra.min()-ra_margin
                      and float(cached.meta.get('RA_MAX',-np.inf)) >= ra.max()+ra_margin
                      and float(cached.meta.get('DEC_MIN',np.inf)) <= dec.min()-margin
                      and float(cached.meta.get('DEC_MAX',-np.inf)) >= dec.max()+margin)
            valid=(int(cached.meta.get('SDSS_DR',-1))==17
                   and set(ns['SDSS_FIELDS']).issubset(cached.colnames)
                   and np.all(np.asarray(cached['mode'])==1)
                   and len(np.unique(cached['objid']))==len(cached))
            if complete and valid:
                print('Reusing SDSS cache: full unchanged 1.5-arcsec match halo is covered.',flush=True)
                return cached
        return original_loader(reference_catalog)

    ns['load_sdss_primary_catalog']=reuse_complete_sdss_cache
    exec(compile(separator.lstrip('\n')+attachment,'singlephot:sdss_attachment','exec'),ns)
    execute(nb, [20,22,24,26,28,30,32,34], ns)
    # Every in-footprint FINALCAT row remains represented in the union.
    baseline = Table.read(BASELINE/f'output/COSMOS{field}/singlephot/cosmos{field}_master_catalog.fits')
    current = ns['master']
    old_idx = np.sort(np.asarray(baseline['finalcat_idx'])[np.asarray(baseline['has_finalcat'],bool)])
    new_idx = np.sort(np.asarray(current['finalcat_idx'])[np.asarray(current['has_finalcat'],bool)])
    assert np.array_equal(old_idx, new_idx), 'FINALCAT footprint changed'
    print(f'COSMOS{field}: complete; FINALCAT footprint unchanged ({len(old_idx):,} rows)', flush=True)


def run_combined():
    nb = notebook('2.3_numofmatchedbands_digonostics.ipynb')
    ns = {'__name__':'combined_flag_ablation'}
    execute(nb, [2], ns)
    output = EXPERIMENT/'singlephot'
    analysis = output/'analysis/numofmatchedbands'
    analysis.mkdir(parents=True, exist_ok=True)
    paths = {field:EXPERIMENT/f'COSMOS{field}/singlephot/cosmos{field}_master_catalog.fits' for field in (1,2,3,4)}
    ns.update(OUTPUT_DIR=output, OVERALL_CATALOG_PATH=output/'master_catalog_MatchedWithFianlcat_OverallCosmosField.fits',
              ANALYSIS_DIR=analysis, EXAMPLE_DIR=analysis/'examples', MASTER_PATHS=paths,
              masters={field:Table.read(path) for field,path in paths.items()})
    # Preserve the original SDSS reconciliation logic and reuse cached queries.
    old_queries = ROOT/'output/singlephot/analysis/numofmatchedbands/sdss_requery'
    if old_queries.exists():
        shutil.copytree(old_queries, analysis/'sdss_requery', dirs_exist_ok=True)
    execute(nb, [4], ns)
    sdss_code=cell_source(nb,6)
    definitions,separator,resolution=sdss_code.partition('# Identify disagreements before resolving any assignments.')
    assert separator
    exec(compile(definitions,'combined:sdss_definitions','exec'),ns)
    remote_query=ns['query_sdss_positions']

    def reuse_local_sdss_candidates(row_indices):
        """Use complete DR17 primary caches for the unchanged candidate query.

        Before replacing a network query, prove every 1.51-arcsec search
        circle lies inside a cached SQL rectangle. Retrieve all cached
        primary objects within those circles, rather than only reusing the
        previous best matches. The downstream disambiguation is unchanged.
        """
        from astropy.coordinates import SkyCoord, search_around_sky
        import astropy.units as u
        targets=ns['finalcat_coordinates'][row_indices]
        radius=(ns['SDSS_MATCH_RADIUS_ARCSEC']+.01)*u.arcsec
        dec=np.asarray(targets.dec.deg)
        ra=np.asarray(targets.ra.deg)
        dec_margin=radius.to_value(u.deg)
        ra_margin=dec_margin/np.cos(np.deg2rad(np.abs(dec)+dec_margin))
        covered=np.zeros(len(row_indices),bool)
        cached_catalogs=[]
        for field in (1,2,3,4):
            cache=EXPERIMENT/f'COSMOS{field}/singlephot/diagnostics/cosmos{field}_master_sdss_matching_primary_dr17.fits'
            table=Table.read(cache)
            assert int(table.meta['SDSS_DR'])==17 and np.all(np.asarray(table['mode'])==1)
            covered|=((ra-ra_margin>=float(table.meta['RA_MIN'])) & (ra+ra_margin<=float(table.meta['RA_MAX']))
                      & (dec-dec_margin>=float(table.meta['DEC_MIN'])) & (dec+dec_margin<=float(table.meta['DEC_MAX'])))
            cached_catalogs.append(table)
        if not covered.all():
            return remote_query(row_indices)
        payloads={}
        for table in cached_catalogs:
            coordinates=SkyCoord(np.asarray(table['ra'],float)*u.deg,np.asarray(table['dec'],float)*u.deg)
            _,selected,_,_=search_around_sky(targets,coordinates,radius)
            for j in np.unique(selected):
                row=table[int(j)]
                payload={name:row[name] for name in ns['SDSS_FIELDS']}
                objid=int(row['objid'])
                if objid in payloads:
                    assert ns['payloads_equal'](payloads[objid],payload), 'SDSS caches disagree'
                payloads[objid]=payload
        print(f'SDSS: resolved candidate retrieval for {len(row_indices)} targets from complete local DR17 caches.',flush=True)
        Table({'master_id':row_indices,'complete_query_circle_in_cache':covered}).write(
            analysis/'sdss_local_candidate_cache_coverage.csv',format='ascii.csv',overwrite=True)
        return payloads

    ns['query_sdss_positions']=reuse_local_sdss_candidates
    exec(compile(separator+resolution,'combined:sdss_resolution','exec'),ns)
    execute(nb, [8,10,12,13,17,19], ns)
    # The cumulative cell is independently runnable and declares its own paths.
    # Redirect only the output and combined-catalog paths; keep FINALCAT input.
    cumulative = cell_source(nb,15)
    # This experiment requests fraction figures; keep optional cumulative
    # output as PNG without generating an additional document artifact.
    cumulative = '\n'.join(line for line in cumulative.splitlines() if "fig.savefig(CUM_OUTPUT / 'rm_cumulative_source_count_by_band.pdf'" not in line)
    cumulative = cumulative.replace("CUM_ROOT / 'output/singlephot/analysis/numofmatchedbands'",repr(analysis))
    cumulative = cumulative.replace("CUM_ROOT / 'output/singlephot/master_catalog_MatchedWithFianlcat_OverallCosmosField.fits'",repr(ns['OVERALL_CATALOG_PATH']))
    # repr(Path) emits PosixPath; inject an explicit constructor for clarity.
    from pathlib import PosixPath
    ns['PosixPath'] = PosixPath
    exec(compile(cumulative,'combined:cumulative','exec'),ns)
    print('Combined catalog and original fraction plots rebuilt.',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['fingerprint','field','combined','verify'])
    parser.add_argument('field',type=int,nargs='?',choices=[1,2,3,4])
    args=parser.parse_args()
    os.chdir(ROOT)
    EXPERIMENT.mkdir(parents=True,exist_ok=True)
    fingerprint_path=EXPERIMENT/'raw_and_mask_fingerprints.json'
    if args.stage=='fingerprint':
        assert not fingerprint_path.exists(), 'Never replace the pre-experiment fingerprints.'
        fingerprint_path.write_text(json.dumps(input_fingerprints(),indent=2)+'\n')
    elif args.stage=='field':
        assert args.field is not None
        run_field(args.field)
    elif args.stage=='combined':
        run_combined()
    else:
        previous=json.loads(fingerprint_path.read_text())
        assert input_fingerprints()==previous, 'Raw catalogs or boundary masks changed.'
        print(f'Verified unchanged hashes, sizes and mtimes for {len(previous)} raw catalogs/masks.')


if __name__=='__main__':
    main()
