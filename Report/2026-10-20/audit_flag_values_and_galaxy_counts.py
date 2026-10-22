"""Audit flag bits and reconstruct galaxy selections from saved masters.

Count newly admitted *single-band measurements* separately from astronomical
objects. A single object has one aperture flag per measured band. Stable raw
source IDs are unique within a field/band only. SDSS objid is used solely for
the comparison of the strict SDSS-galaxy subsets, never FINALCAT ID strings.
"""
from pathlib import Path
import warnings
from collections import Counter
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.units import UnitsWarning
import nbformat

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / 'output/experiments/aperture_flag_ablation'
OUT = EXP / 'statistical_analysis'


def val(t, name, fill=np.nan):
    return np.asarray(np.ma.filled(t[name], fill))


def strict_galaxies(t, bands):
    """Reproduce notebook S3, starting from the complete master itself."""
    has = val(t, 'has_7dt', False).astype(bool)
    galaxy = has & val(t, 'has_sdss', False).astype(bool) & (val(t, 'sdss_type', -1) == 3)
    phot = galaxy.copy()
    for band in 'gri':
        mag = val(t, f'sdss_petroMag_{band}')
        phot &= np.isfinite(mag) & (mag > -50) & (mag < 50)
    clean = phot & (val(t, 'sdss_clean', 0) == 1) & (val(t, 'sdss_mode', 0) == 1)
    unique = clean & (val(t, 'sdss_n_candidates', 0) == 1) & (val(t, 'sdss_n_master_candidates', 0) == 1)
    full = unique & (val(t, 'n_7dt_bands', 0) == len(bands))
    for band in bands:
        mag, err = val(t, band+'_mag'), val(t, band+'_mag_err')
        full &= np.isfinite(mag) & (mag > 5) & np.isfinite(err) & (err > 0)
    return dict(any_7dt_sdss_galaxy=galaxy, valid_sdss_gri=phot,
                clean_primary=clean, clean_primary_one_to_one=unique,
                full_available_bands=full,
                any_7dt_sdss_galaxy_also_finalcat=galaxy & val(t,'has_finalcat',False).astype(bool))


def main():
    warnings.filterwarnings('ignore', category=UnitsWarning)
    OUT.mkdir(parents=True, exist_ok=True)
    nb = nbformat.read(ROOT/'Notebooks/2_singphotandmatching.ipynb',4)
    # Reuse the actual selection function; provide only its mask dependency.
    # This executes definitions only and does not load science images or run SE++.
    source = nb.cells[2].source
    definitions = source[source.index('def positions_in_valid_mask'):source.index('\n\nprint(f"Project root:')]
    selected_rows, stage_rows, signature_rows, galaxy_flag_rows = [], [], [], []
    pooled = Counter()
    for f in range(1,5):
        d = ROOT/f'output/COSMOS{f}/singlephot'
        masters = [Table.read(EXP/f'baseline/output/COSMOS{f}/singlephot/cosmos{f}_master_catalog.fits'),
                   Table.read(d/f'cosmos{f}_master_catalog.fits')]
        bands = masters[0].meta['ACTBANDS'].split(',')
        namespace = dict(np=np, valid_mask=fits.getdata(ROOT/f'COSMOS_data/COSMOS{f}/COSMOS{f}_validmask.fits').astype(bool))
        exec(compile(definitions,'selection definitions','exec'),namespace)
        masks = [strict_galaxies(t,bands) for t in masters]
        for label,t,selection in zip(['before','after'],masters,masks):
            for stage,use in selection.items():
                stage_rows.append(dict(COSMOS=f,variant=label,stage=stage,N_master_rows=int(use.sum()),
                    N_unique_SDSS_objid=len(np.unique(val(t,'sdss_objid',-1)[use]))))
            # Independently reconstructed S3 must agree with the saved diagnostic.
            diag = (EXP/f'baseline/full_products/COSMOS{f}/singlephot/diagnostics') if label=='before' else d/'diagnostics'
            summary = Table.read(diag/f'cosmos{f}_7dt_positional_quality_sample_summary.fits')
            stored = summary[np.asarray(summary['sample']).astype(str)=='S3_galaxy']['n_sources'][0]
            assert int(selection['full_available_bands'].sum()) == int(stored)
        totals = Counter()
        flags_by_source = {}
        for band in bands:
            raw = Table.read(d/f'cat_cosmos{f}_{band}.fits')
            chosen = namespace['select_singleband_sources'](raw)
            flags = val(chosen,'aper3_flags',-1).astype(int)
            counter = Counter(flags[flags!=0].tolist())
            totals.update(counter)
            for flag,n in sorted(counter.items()):
                selected_rows.append(dict(COSMOS=f,band=band,flag=flag,N_added_measurements=n))
            flags_by_source[band] = {int(i):int(flag) for i,flag in zip(val(raw,'source_id',-1),val(raw,'aper3_flags',-1))}
        selection_check = Table.read(EXP/f'COSMOS{f}/singlephot/selection_comparison.csv')
        assert sum(totals.values()) == int(np.sum(selection_check['n_added']))
        pooled.update(totals)
        for flag,n in sorted(totals.items()):
            selected_rows.append(dict(COSMOS=f,band='ALL',flag=flag,N_added_measurements=n))
        if f == 1:
            old_ids = set(val(masters[0],'sdss_objid',-1)[masks[0]['full_available_bands']].tolist())
            new_ids = set(val(masters[1],'sdss_objid',-1)[masks[1]['full_available_bands']].tolist())
            after = masters[1]
            full = val(after,'has_7dt',False).astype(bool) & (val(after,'n_7dt_bands',0)==20)
            for band in bands:
                mag,err=val(after,band+'_mag'),val(after,band+'_mag_err')
                full &= np.isfinite(mag)&(mag>5)&np.isfinite(err)&(err>0)
            old_anchor = set(val(masters[0],'m600_source_id',-1)[val(masters[0],'n_7dt_bands',0)==20].tolist())
            cases = [
                ('new_full20_sources', full & ~np.isin(val(after,'m600_source_id',-1),list(old_anchor))),
                ('new_strict_SDSS_full20_galaxies', masks[1]['full_available_bands'] & ~np.isin(val(after,'sdss_objid',-1),list(old_ids)))]
            for name,use in cases:
                measures = Counter(); signatures = Counter(); nflagged=[]
                for i in np.flatnonzero(use):
                    perband = [flags_by_source[b][int(after[b+'_source_id'][i])] for b in bands]
                    measures.update(x for x in perband if x!=0)
                    signatures[tuple(sorted(set(perband)-{0}))]+=1
                    nflagged.append(sum(x!=0 for x in perband))
                    if name=='new_strict_SDSS_full20_galaxies':
                        galaxy_flag_rows.append(dict(sdss_objid=int(after['sdss_objid'][i]),
                            branch_master_id=int(after['master_id'][i]),
                            N_nonzero_flag_bands=nflagged[-1],
                            **{b+'_flag':flag for b,flag in zip(bands,perband)}))
                for flag,n in sorted(measures.items()):
                    signature_rows.append(dict(population=name,kind='flagged_band_measurements',flag_or_signature=str(flag),N=n))
                for sig,n in sorted(signatures.items()):
                    signature_rows.append(dict(population=name,kind='objects_by_flag_values',flag_or_signature='+'.join(map(str,sig)),N=n))
                signature_rows.append(dict(population=name,kind='objects_total',flag_or_signature='all',N=int(use.sum())))
                print(name, 'N',int(use.sum()), 'flagged-band median/range',np.median(nflagged),min(nflagged),max(nflagged),'flags',dict(measures),'signatures',dict(signatures))
            print('Strict galaxies old/new/common/gained/lost:',len(old_ids),len(new_ids),len(old_ids&new_ids),len(new_ids-old_ids),len(old_ids-new_ids))
    for flag,n in sorted(pooled.items()):
        selected_rows.append(dict(COSMOS=0,band='ALL',flag=flag,N_added_measurements=n))
    tables = dict(added_measurement_flag_distribution=Table(rows=selected_rows),
                  sdss_galaxy_selection_stages=Table(rows=stage_rows),
                  cosmos1_new_fullband_flag_populations=Table(rows=signature_rows),
                  cosmos1_new_strict_full20_galaxy_flags=Table(rows=galaxy_flag_rows))
    for name,t in tables.items():
        t.write(OUT/(name+'.csv'),format='ascii.csv',overwrite=True)
    print('Added measurement flags by field (COSMOS=0: all fields)')
    t=tables['added_measurement_flag_distribution']; t[t['band']=='ALL'].pprint(max_lines=50,max_width=180)
    print('Galaxy stages COSMOS1: derived directly from the two master catalogs')
    t=tables['sdss_galaxy_selection_stages'];t[t['COSMOS']==1].pprint(max_lines=30,max_width=180)
    return OUT


def save_notebook_outputs():
    """Preserve the selection clarification beside the preceding audit."""
    from nbclient import NotebookClient
    heading = '## Flag meanings and clarification of the 10-to-63 galaxy subset'
    prose = heading + '''

The 10 baseline galaxies are reconstructed from this experiment's saved
pre-ablation master (`d3c9272` baseline); this count was not copied from the
historical Proposal. It is **not the number of all SDSS galaxy matches**.
The selection-stage table distinguishes any-band SDSS galaxies, valid gri,
clean/primary matches, one-to-one matches, and complete available-band S3.

Installed SE++ `SEFramework/Source/SourceFlags.h` agrees with the
[official flag enum](https://github.com/astrorama/SourceXtractorPlusPlus/blob/develop/SEFramework/SEFramework/Source/SourceFlags.h):
0 = no flags set; 2 = BLENDED; 4 = SATURATED; 6 = BLENDED | SATURATED;
16 = NEIGHBORS; 18 = BLENDED | NEIGHBORS. These are bit masks, not ordered
quality scores. In particular, SE++ bit 16 must not be interpreted using the
legacy SExtractor flag table. Nonzero flags identify a condition to inspect,
not a universal declaration that the flux is unusable; zero is not an accuracy
guarantee. Removing the complete flag cut also admits saturation warnings.

The measurement histogram counts newly admitted raw **field-band detections**
after all unchanged finite-value, radius, magnitude and footprint cuts.
It does not count unique objects. Each object has a separate flag in each band.
The COSMOS1 galaxy comparison uses SDSS objid to distinguish retained/new
members of the strict complete-20-band S3 subset. All 10 baseline members
remain, with 53 additional members. Object-signature rows list the set of
nonzero flag values encountered across that object's bands.
'''
    code = '''# Read-only audit of saved raw catalogs and both master reductions.
from pathlib import Path
import importlib.util
from astropy.table import Table
from IPython.display import display
FLAG_ROOT = Path.cwd().resolve()
if FLAG_ROOT.name == 'Notebooks':
    FLAG_ROOT = FLAG_ROOT.parent
FLAG_SPEC = importlib.util.spec_from_file_location(
    'flag_values_and_galaxy_counts',
    FLAG_ROOT / 'Report/2026-10-20/audit_flag_values_and_galaxy_counts.py')
FLAG_MODULE = importlib.util.module_from_spec(FLAG_SPEC)
FLAG_SPEC.loader.exec_module(FLAG_MODULE)
FLAG_OUT = FLAG_MODULE.main()
for FLAG_NAME in ['added_measurement_flag_distribution',
                  'sdss_galaxy_selection_stages',
                  'cosmos1_new_fullband_flag_populations']:
    FLAG_TABLE = Table.read(FLAG_OUT / (FLAG_NAME + '.csv'))
    if FLAG_NAME == 'added_measurement_flag_distribution':
        FLAG_TABLE = FLAG_TABLE[FLAG_TABLE['band'] == 'ALL']
    print(FLAG_NAME)
    display(FLAG_TABLE)
'''
    audit = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(code)])
    audit.metadata.kernelspec = dict(name='python3', display_name='Python 3', language='python')
    executed = NotebookClient(audit, timeout=600).execute(cwd=str(ROOT))
    for name in ['Proposal.ipynb', '4_Proposal_aperture_flag_results.ipynb',
                 '2.3_numofmatchedbands_digonostics.ipynb']:
        path = ROOT/'Notebooks'/name
        notebook = nbformat.read(path,4)
        for i,cell in enumerate(notebook.cells):
            if cell.cell_type=='markdown' and cell.source.startswith(heading):
                del notebook.cells[i:i+2]
                break
        notebook.cells.extend([nbformat.v4.new_markdown_cell(prose),executed.cells[0]])
        nbformat.validate(notebook)
        nbformat.write(notebook,path)
        print('Saved flag/selection clarification:',path)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--update-notebooks', action='store_true')
    args = parser.parse_args()
    if args.update_notebooks:
        save_notebook_outputs()
    else:
        main()
