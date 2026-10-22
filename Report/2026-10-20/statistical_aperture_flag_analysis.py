"""Paired catalog audit of the aperture-flag ablation; never run photometry.

FINALCAT row indices identify reference objects. Local master IDs and FINALCAT
ID strings must not be joined across reductions: the former are regenerated
and the latter are not unique. Positional residuals are summarized per object
before comparing populations, so 20 bands are not treated as 20 independent
galaxies. The random-offset diagnostic is a chance-alignment proxy, not a
calibrated false-discovery rate or a test of photometric-redshift accuracy.
"""
from pathlib import Path
import importlib.util
import json
import warnings
import numpy as np
from astropy.table import Table
from astropy.units import UnitsWarning
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / 'output/experiments/aperture_flag_ablation'
OUT = EXP / 'statistical_analysis'
BANDS = tuple(f'm{w}' for w in range(400, 900, 25))


def values(table, column, fill=np.nan):
    """Respect FITS masks rather than interpreting masked backing values."""
    return np.asarray(np.ma.filled(table[column], fill))


def save(rows, name):
    table = Table(rows=rows)
    table.write(OUT / f'{name}.csv', format='ascii.csv', overwrite=True)
    return table


def main():
    # Raw catalogs carry a legacy nonstandard unit string. Numeric values used
    # here do not depend on that unit, so suppress only this parser warning.
    warnings.filterwarnings('ignore', category=UnitsWarning)
    OUT.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('fraction_helper', Path(__file__).with_name('plot_aperture_flag_comparison.py'))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    filename = helper.FILENAME
    reference = Table.read(ROOT / 'COSMOS_data/COSMOS_finalcat.fits')
    rm = values(reference, 'RM')
    valid = np.isfinite(rm) & (rm > 0) & (rm < 50)
    bright = valid & (rm <= 21)
    samples = [helper.selected_columns(EXP / 'baseline/output/singlephot' / filename),
               helper.selected_columns(ROOT / 'output/singlephot' / filename)]
    counts = []
    for sample in samples:
        assert np.all((sample['ids'] >= 0) & (sample['ids'] < len(reference)))
        assert np.allclose(sample['rm'], rm[sample['ids']], equal_nan=True)
        n = np.zeros(len(reference), dtype=int)
        n[sample['ids']] = sample['nband']
        counts.append(n)
    before, after = counts
    delta = after - before
    coverage = Table.read(ROOT / 'output/singlephot/diagnostics/cosmos1-4_7dt_finalcat_matching_coverage_catalog.fits')
    available = np.zeros((len(reference), 20), dtype=bool)
    for f in range(1, 5):
        meta = Table.read(ROOT / f'output/COSMOS{f}/singlephot/cosmos{f}_master_catalog.fits').meta
        active = np.array([b in meta['ACTBANDS'].split(',') for b in BANDS])
        # Coverage is source-level footprint eligibility, not a guarantee that
        # each band has usable pixels or sufficient depth at that position.
        available[values(coverage, 'finalcat_idx', -1)] |= values(coverage, f'COSMOS{f}_in_coverage', False)[:, None] & active
    n_available = available.sum(axis=1)
    assert np.all(n_available[after == 20] == 20)
    populations = []
    for label, use in [('all_RMle21', bright), ('footprint_RMle21', bright & (n_available > 0)),
                       ('20band_eligible_RMle21', bright & (n_available == 20))]:
        for name, threshold in [('any', 1), ('ge15', 15), ('all20', 20)]:
            a, b = before[use] >= threshold, after[use] >= threshold
            populations.append(dict(sample=label, population=name, N=int(use.sum()),
                before=int(a.sum()), after=int(b.sum()), gained=int((~a & b).sum()),
                lost=int((a & ~b).sum()), before_percent=100*a.mean(),
                after_percent=100*b.mean(), difference_pp=100*(b.mean()-a.mean())))
    pop = save(populations, 'paired_recovery')
    strata = []
    for lo, hi in [(0,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,50)]:
        use = valid & (rm > lo) & (rm <= hi)
        row = dict(RM_lo=lo, RM_hi=hi, N=int(use.sum()), median_bands_before=float(np.median(before[use])),
                   median_bands_after=float(np.median(after[use])), gained_bands=int(delta[use].sum()),
                   improved=int((delta[use]>0).sum()), unchanged=int((delta[use]==0).sum()), worsened=int((delta[use]<0).sum()))
        for label, n in [('before',before),('after',after)]:
            for name, threshold in [('any',1),('ge15',15),('all20',20)]:
                row[f'{name}_{label}_percent'] = float(100*np.mean(n[use]>=threshold))
        strata.append(row)
    magnitude = save(strata, 'magnitude_strata')
    save([dict(sample=name, N=int(use.sum()), improved=int((delta[use]>0).sum()),
               unchanged=int((delta[use]==0).sum()), worsened=int((delta[use]<0).sum()),
               median_delta_bands=float(np.median(delta[use])), total_delta_bands=int(delta[use].sum()))
          for name,use in [('RMle21',bright),('RMle21_baseline_detected',bright&(before>0)),
                           ('RMle21_detected_in_either',bright&((before>0)|(after>0)))]], 'paired_band_changes')

    positional, radius, raw_quality, assignments, selected_types = [], [], [], [], []
    for f in range(1,5):
        directories = [EXP/f'baseline/output/COSMOS{f}/singlephot', ROOT/f'output/COSMOS{f}/singlephot']
        diagnostics = [EXP/f'baseline/full_products/COSMOS{f}/singlephot/diagnostics', directories[1]/'diagnostics']
        masters = [Table.read(d/f'cosmos{f}_master_catalog.fits') for d in directories]
        matched = []
        for t in masters:
            use = values(t, 'has_7dt', False) & values(t, 'has_finalcat', False)
            matched.append({int(t['finalcat_idx'][i]):i for i in np.flatnonzero(use)})
        common = sorted(set(matched[0]) & set(matched[1]))
        same, changed, err_delta, mag_delta = 0, 0, [], []
        for band in masters[0].meta['ACTBANDS'].split(','):
            ia, ib = ([m[k] for k in common] for m in matched)
            sa = values(masters[0], band+'_source_id', -1)[ia]
            sb = values(masters[1], band+'_source_id', -1)[ib]
            both = (sa>=0)&(sb>=0)
            unchanged = both & (sa==sb)
            same += int(unchanged.sum()); changed += int((both & (sa!=sb)).sum())
            for suffix, target in [('_mag',mag_delta),('_mag_err',err_delta)]:
                a = values(masters[0],band+suffix)[ia]
                b = values(masters[1],band+suffix)[ib]
                target.extend((b-a)[unchanged].tolist())
            raw = Table.read(directories[1]/f'cat_cosmos{f}_{band}.fits')
            order = np.argsort(values(raw,'source_id',-1))
            ids = values(raw,'source_id',-1)[order]
            source_sets = [set(values(t,band+'_source_id',-1)[values(t,band+'_source_id',-1)>=0]) for t in masters]
            for name, selected in [('retained',source_sets[0]&source_sets[1]), ('added',source_sets[1]-source_sets[0])]:
                rows = order[np.searchsorted(ids,sorted(selected))]
                assert np.array_equal(values(raw,'source_id',-1)[rows], sorted(selected))
                errors = values(raw,'aper3_mag_err')[rows]
                snr = values(raw,'snrratio')[rows]
                flags = values(raw,'aper3_flags',-1)[rows]
                raw_quality.append(dict(COSMOS=f,band=band,population=name,N=len(rows),
                    median_mag_err=float(np.nanmedian(errors)), median_detection_snr=float(np.nanmedian(snr)),
                    frac_mag_err_gt_02=float(np.mean(errors>.2)), frac_nonzero_flag=float(np.mean(flags!=0))))
        assignments.append(dict(COSMOS=f,common_reference_sources=len(common), shared_band_associations=same,
            changed_band_associations=changed, fraction_changed_among_both=changed/(same+changed),
            max_abs_mag_change_same_raw=float(np.nanmax(np.abs(mag_delta))),
            max_abs_magerr_change_same_raw=float(np.nanmax(np.abs(err_delta)))))

        # Use the stable m600 raw detection as an object anchor across runs.
        # Require every active-band source ID to agree for the paired control;
        # changed memberships are reported separately and never called controls.
        grouped = []
        for diag in diagnostics:
            t = Table.read(diag/f'cosmos{f}_7dt_positional_quality_fullband_catalog.fits')
            groups = {}
            for mid in np.unique(t['master_id']):
                sub = t[t['master_id']==mid]
                anchor = sub[np.asarray(sub['band']).astype(str)=='m600']
                if len(anchor)!=1: continue
                anchor_id = int(anchor['source_id'][0])
                assert anchor_id not in groups
                groups[anchor_id] = dict(signature=tuple(sorted(zip(sub['band'],sub['source_id']))),
                    residual=float(np.nanmedian(sub['loo_residual_arcsec'])),
                    snr=float(np.nanmedian(sub['snrratio'])))
            grouped.append(groups)
        old, new = grouped
        common_anchor = set(old)&set(new)
        stable = [k for k in common_anchor if old[k]['signature']==new[k]['signature']]
        for name, keys, source in [('baseline_all',list(old),old),('branch_all',list(new),new),
                                   ('stable_before',stable,old),('stable_after',stable,new),
                                   ('new_fullband',list(set(new)-set(old)),new)]:
            positional.append(dict(COSMOS=f,population=name,N=len(keys),
                median_object_loo_arcsec=float(np.median([source[k]['residual'] for k in keys])),
                median_object_snr=float(np.median([source[k]['snr'] for k in keys])),
                stable_pairs=len(stable), changed_memberships=len(common_anchor)-len(stable)))
        for label,diag in zip(['before','after'],diagnostics):
            t = Table.read(diag/f'cosmos{f}_7dt_positional_quality_sample_summary.fits')
            for row in t:
                selected_types.append(dict(COSMOS=f,variant=label,sample=str(row['sample']),N=int(row['n_sources'])))
            t = Table.read(diag/f'cosmos{f}_7dt_finalcat_matching_radius_diagnostics.fits')
            t = t[np.isclose(t['radius_arcsec'],.7)]
            actual = t[np.asarray(t['sample']).astype(str)=='actual'][0]
            random = np.asarray(t['n_selected'][np.asarray(t['sample']).astype(str)=='random'])
            radius.append(dict(COSMOS=f,variant=label,N_actual=int(actual['n_selected']),
                median_sep_arcsec=float(actual['median_sep_arcsec']), random_mean=float(random.mean()),
                random_min=int(random.min()),random_max=int(random.max()),
                random_over_actual_percent=100*float(random.mean()/actual['n_selected'])))
    position = save(positional,'positional_population_control')
    radii = save(radius,'chance_alignment_proxy')
    save(raw_quality,'raw_measurement_quality'); save(assignments,'same_object_associations')
    save(selected_types,'fullband_sample_composition')
    depth = Table.read(EXP/'comparison/depth_before_after.csv')
    depth_delta = values(depth,'branch_median_depth_AB')-values(depth,'baseline_median_depth_AB')
    summary = dict(photometry_rerun=False, N_reference_bright=int(bright.sum()),
        N_bright_20band_eligible=int((bright&(n_available==20)).sum()),
        depth_median_shift_mag=float(np.nanmedian(depth_delta)),
        depth_shift_min_mag=float(np.nanmin(depth_delta)),depth_shift_max_mag=float(np.nanmax(depth_delta)),
        explanation='Depth changes reflect changed catalog-derived aperture size, not new image information. Recovery and random offsets do not measure photo-z accuracy.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

    fig, axes = plt.subplots(2,2,figsize=(12,9),layout='constrained')
    x = np.arange(len(magnitude)-1)
    labels = ['<18','18–19','19–20','20–21','21–22','22–23']
    for name,color in [('any','C0'),('ge15','C1'),('all20','C2')]:
        for label,style in [('before','--'),('after','-')]:
            axes[0,0].plot(x,magnitude[f'{name}_{label}_percent'][:-1],style+'o',color=color,label=f'{name}: {label}')
    axes[0,0].set(ylabel='Full-reference recovery [%]',xlabel='FINALCAT RM bin')
    axes[0,0].set_xticks(x,labels); axes[0,0].legend(fontsize=8,ncol=2)
    use = bright & ((before>0)|(after>0))
    axes[0,1].hist(delta[use],bins=np.arange(-20.5,21),color='C0')
    axes[0,1].set(xlabel='Change in successful-band count, same reference object',ylabel='Sources (RM <= 21; detected in either)')
    axes[0,1].axvline(0,color='k',lw=1)
    for name,offset,color in [('stable_before',-.2,'C0'),('stable_after',0,'C1'),('new_fullband',.2,'C2')]:
        sub=position[position['population']==name]
        axes[1,0].bar(np.arange(1,5)+offset,sub['median_object_loo_arcsec'],width=.2,label=name,color=color)
    axes[1,0].set(xlabel='COSMOS field',ylabel='Median per-object LOO residual [arcsec]'); axes[1,0].set_xticks(range(1,5));axes[1,0].legend(fontsize=8)
    for label,offset,color in [('before',-.15,'C0'),('after',.15,'C1')]:
        sub=radii[radii['variant']==label]
        axes[1,1].bar(np.arange(1,5)+offset,sub['random_over_actual_percent'],width=.3,label=label,color=color)
    axes[1,1].set(xlabel='COSMOS field',ylabel='Offset-match count / actual [%] (proxy)');axes[1,1].set_xticks(range(1,5));axes[1,1].legend()
    for ax in axes.flat: ax.grid(axis='y',alpha=.2)
    fig.suptitle('Aperture-flag ablation: recovery, paired controls and chance-alignment proxy')
    fig.savefig(OUT/'statistical_audit.png',dpi=170)
    plt.close(fig)
    print(pop); print(magnitude); print(position); print(radii); print(summary)
    return OUT


def save_notebook_outputs():
    """Execute a self-contained audit cell and preserve real rendered outputs.

    The audit reads existing products only. The isolated kernel avoids rerunning
    upstream notebook cells that could invoke image-level analysis.
    """
    import nbformat
    from nbclient import NotebookClient
    heading = '## Statistical interpretation of the aperture-flag ablation'
    prose = heading + '''

This is a paired comparison of the same FINALCAT rows, not two independent
samples. Report gained and lost objects separately. Fractions use the unchanged
full-reference denominator; the eligible-footprint rows additionally show the
effect of missing bands. Eligibility means field coverage with an active band,
not guaranteed usable pixels or adequate depth at each position.

For positional quality, aggregate leave-one-band-out (LOO) residuals per object
first. The paired control requires exactly the same raw detection in every
active band, anchored by the stable m600 detection ID. Differences between the
overall medians can reflect population composition rather than changed raw
photometry. S1/S3 are clean, primary, one-to-one SDSS stars/galaxies with valid
gri and all available 7DT bands; they are selected subsets, not complete class
censuses. Offset matches at 0.7 arcsec are a chance-alignment proxy, not measured
contamination or photo-z accuracy. These deterministic catalog differences are
effect sizes; an independent-sample significance test would not answer whether
the additional photometry is scientifically accurate.
'''
    code = '''# Audit existing catalogs only; no photometry or photo-z fitting.
from pathlib import Path
import importlib.util
from astropy.table import Table
from IPython.display import display, Image

AUDIT_ROOT = Path.cwd().resolve()
if AUDIT_ROOT.name == 'Notebooks':
    AUDIT_ROOT = AUDIT_ROOT.parent
AUDIT_SPEC = importlib.util.spec_from_file_location(
    'flag_statistical_audit',
    AUDIT_ROOT / 'Report/2026-10-20/statistical_aperture_flag_analysis.py')
AUDIT_MODULE = importlib.util.module_from_spec(AUDIT_SPEC)
AUDIT_SPEC.loader.exec_module(AUDIT_MODULE)
AUDIT_OUT = AUDIT_MODULE.main()
display(Image(filename=str(AUDIT_OUT / 'statistical_audit.png')))
for AUDIT_NAME in ['paired_recovery', 'paired_band_changes', 'magnitude_strata',
                   'positional_population_control', 'same_object_associations',
                   'fullband_sample_composition', 'chance_alignment_proxy']:
    print(AUDIT_NAME)
    display(Table.read(AUDIT_OUT / (AUDIT_NAME + '.csv')))
'''
    audit = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(code)])
    audit.metadata.kernelspec = dict(name='python3', display_name='Python 3', language='python')
    executed = NotebookClient(audit, timeout=600).execute(cwd=str(ROOT))
    for name in ['Proposal.ipynb','4_Proposal_aperture_flag_results.ipynb',
                 '2.3_numofmatchedbands_digonostics.ipynb']:
        path = ROOT / 'Notebooks' / name
        notebook = nbformat.read(path, as_version=4)
        for i, cell in enumerate(notebook.cells):
            if cell.cell_type == 'markdown' and cell.source.startswith(heading):
                del notebook.cells[i:i+2]
                break
        notebook.cells.extend([nbformat.v4.new_markdown_cell(prose), executed.cells[0]])
        nbformat.validate(notebook)
        nbformat.write(notebook, path)
        print(f'Saved executed statistical audit: {path}')


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--update-notebooks', action='store_true')
    arguments = parser.parse_args()
    if arguments.update_notebooks:
        save_notebook_outputs()
    else:
        main()
