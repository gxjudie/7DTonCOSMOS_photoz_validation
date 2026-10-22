"""Export two clearly defined subsets per COSMOS field directly under output/.

Local master IDs are regenerated, so identify pre-existing 7DT objects by
shared raw source IDs within the SAME band. FINALCAT row indices independently
identify reference objects; FINALCAT string IDs need not be unique. An old
FINALCAT-only row is not evidence that a 7DT detection existed previously.
Full-band means all available field bands: COSMOS1=20, 2=16, 3=18, 4=16.
"""
from collections import Counter
from pathlib import Path
import warnings
import numpy as np
from astropy.table import Table, MaskedColumn
from astropy.units import UnitsWarning

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT/'output/experiments/aperture_flag_ablation'
BANDS = tuple(f'm{w}' for w in range(400,900,25))


def values(table, name, fill=-1):
    return np.asarray(np.ma.filled(table[name], fill))


def export_field(field):
    """Retain the same object-identity and selection rules in every field."""
    warnings.filterwarnings('ignore',category=UnitsWarning)
    current_path = ROOT/f'output/COSMOS{field}/singlephot/cosmos{field}_master_catalog.fits'
    baseline_path = EXP/f'baseline/output/COSMOS{field}/singlephot/cosmos{field}_master_catalog.fits'
    before, after = Table.read(baseline_path), Table.read(current_path)
    active_bands = tuple(after.meta['ACTBANDS'].split(','))
    assert active_bands == tuple(before.meta['ACTBANDS'].split(','))
    n_active = len(active_bands)
    # Some field catalogs omit entirely unobserved-band columns. Keep nominal
    # 20-band audit columns, with absent measurements masked, in every export.
    def raw_ids(table):
        return np.column_stack([values(table,b+'_source_id') if b+'_source_id' in table.colnames
                                else np.full(len(table),-1,dtype=np.int64) for b in BANDS])
    raw_before, raw_after = raw_ids(before), raw_ids(after)
    lookup = [{int(sid):int(i) for i,sid in enumerate(raw_before[:,j]) if sid>=0}
              for j in range(20)]
    old_counts = values(before,'n_7dt_bands',0)
    old_ids = values(before,'master_id')
    ref_lookup = {int(before['finalcat_idx'][i]):int(i)
                  for i in np.flatnonzero(np.asarray(before['has_finalcat'],bool))}
    current_full = np.asarray(after['has_7dt'],bool) & (values(after,'n_7dt_bands',0)==n_active)
    old_full = np.asarray(before['has_7dt'],bool) & (old_counts==n_active)
    # Match the existing per-field complete-band diagnostic definition.
    for band in active_bands:
        for table,selection in [(before,old_full),(after,current_full)]:
            mag = values(table,band+'_mag',np.nan)
            err = values(table,band+'_mag_err',np.nan)
            selection &= np.isfinite(mag)&(mag>5)&np.isfinite(err)&(err>0)
    matched = np.asarray(after['master_class']).astype(str)=='matched'
    candidate_indices = np.flatnonzero(current_full|matched)
    records = {}
    new_full, new_matched = [], []
    for i in candidate_indices:
        support = Counter(lookup[j][int(sid)] for j,sid in enumerate(raw_after[i])
                          if sid>=0 and int(sid) in lookup[j])
        ranked = sorted(support,key=lambda k:(-support[k],-old_counts[k],old_ids[k]))
        parent = ranked[0] if ranked else None
        finalcat_idx = after['finalcat_idx'][i]
        ref_parent = None if np.ma.is_masked(finalcat_idx) else ref_lookup.get(int(finalcat_idx))
        records[int(i)] = (parent,ref_parent,ranked,support)
        # Include objects with no prior 7DT detections (0->all available).
        # Exclude any source linked to an already-complete old object.
        if current_full[i] and not any(old_full[k] for k in ranked):
            new_full.append(int(i))
        # Strictly new matched detections: no retained raw measurement in any
        # band AND the reference row previously had no 7DT counterpart. This
        # excludes old 7DT-only sources that merely acquire a FINALCAT match,
        # and replacement detections for an already matched reference object.
        if matched[i] and not support and ref_parent is not None and not bool(before['has_7dt'][ref_parent]):
            new_matched.append(int(i))

    raw_flag_lookup = []
    for band in BANDS:
        if band in active_bands:
            raw = Table.read(ROOT/f'output/COSMOS{field}/singlephot/cat_cosmos{field}_{band}.fits')
            raw_flag_lookup.append(dict(zip(values(raw,'source_id'),values(raw,'aper3_flags'))))
        else:
            raw_flag_lookup.append({})

    full_label = '20band' if n_active==20 else f'full{n_active}band'
    names = (f'COSMOS{field}_flagcut_removed_newly_{full_label}_sources.fits',
             f'COSMOS{field}_flagcut_removed_new_7DT_sources_matched_to_FINALCAT.fits')
    outputs = []
    for selected,name,kind in [(new_full,names[0],'NEW_FULL_AVAILABLE'),(new_matched,names[1],'NEW_MATCHED_7DT')]:
        result = after[selected].copy()
        parents = [records[i][0] for i in selected]
        ref_parents = [records[i][1] for i in selected]
        missing_parent = np.array([p is None for p in parents])
        missing_reference = np.array([p is None for p in ref_parents])
        result['baseline_master_id'] = MaskedColumn([int(old_ids[p]) if p is not None else -1 for p in parents],mask=missing_parent)
        result['baseline_master_id'].description = f'Old 7DT object with most identical raw band detections; IDs are local to COSMOS{field}.'
        result['baseline_n_7dt_bands'] = [int(old_counts[p]) if p is not None else 0 for p in parents]
        result['baseline_n_7dt_bands'].description = 'Band count of selected old 7DT parent; zero means no shared old raw detection.'
        result['baseline_n_source_candidates'] = [len(records[i][2]) for i in selected]
        result['baseline_shared_raw_bands'] = [int(records[i][3].get(p,0)) for i,p in zip(selected,parents)]
        result['baseline_candidate_master_ids'] = MaskedColumn(
            [','.join(str(int(old_ids[k])) for k in records[i][2]) for i in selected],mask=missing_parent)
        result['baseline_finalcat_master_id'] = MaskedColumn([int(old_ids[p]) if p is not None else -1 for p in ref_parents],mask=missing_reference)
        result['baseline_finalcat_master_id'].description = 'Old master row for the same unique FINALCAT row index, including FINALCAT-only rows.'
        result['baseline_finalcat_had_7dt'] = [bool(before['has_7dt'][p]) if p is not None else False for p in ref_parents]
        result['is_entirely_new_7dt_source'] = missing_parent
        lists = [[] for _ in selected]
        new_measurement_count = np.zeros(len(selected),dtype=int)
        for j,band in enumerate(BANDS):
            old_source = np.array([int(raw_before[p,j]) if p is not None else -1 for p in parents])
            current_source = raw_after[selected,j]
            added = (current_source>=0)&(old_source<0)
            is_new_raw = np.array([sid>=0 and int(sid) not in lookup[j] for sid in current_source])
            flag = np.array([int(raw_flag_lookup[j][sid]) if sid>=0 else -1 for sid in current_source],dtype=np.int64)
            result['baseline_'+band+'_source_id'] = MaskedColumn(old_source,mask=old_source<0)
            result[band+'_added_band'] = added
            result[band+'_added_band'].description = 'Current detection fills a band absent from the chosen old 7DT parent.'
            result[band+'_new_raw_detection'] = is_new_raw
            result[band+'_new_raw_detection'].description = 'This exact raw detection did not occur in any old 7DT master row in this band.'
            result[band+'_aper3_flags'] = MaskedColumn(flag,mask=current_source<0)
            result[band+'_aper3_flags'].description = 'Raw SE++ 3-arcsec aperture flag bitmask: 2 blend, 4 saturated, 16 neighbors; sums combine bits.'
            new_measurement_count += is_new_raw
            for k in np.flatnonzero(added): lists[k].append(band)
        result['n_added_bands'] = [len(bands) for bands in lists]
        result['added_bands'] = [','.join(bands) for bands in lists]
        result['n_new_raw_detections'] = new_measurement_count
        result.meta.update(PRODUCT=kind,BASECOM='d3c9272',PARENT=str(current_path.relative_to(ROOT)),
                           BASEFILE=str(baseline_path.relative_to(ROOT)),FLAGCUT='REMOVED',
                           PHOTRUN=False,ID_SCOPE=f'COSMOS{field} field-local master_id',
                           N_ACTIVE=n_active,ACTBANDS=','.join(active_bands),
                           JOINKEY='Raw source_id within same band; independent FINALCAT row index',
                           SELECT=(f'{n_active} valid bands now; no shared old complete-{n_active}-band parent. Includes 0->{n_active}.'
                                   if kind=='NEW_FULL_AVAILABLE' else
                                   'matched now; zero shared old raw detections; same FINALCAT row was FINALCAT-only.'),
                           S3_ONLY=False,NR_SAVED=len(result))
        path = ROOT/'output'/name
        result.write(path,format='fits',overwrite=True)
        restored = Table.read(path)
        assert len(restored)==len(result)
        for column in result.colnames:
            assert np.array_equal(np.ma.getmaskarray(result[column]),np.ma.getmaskarray(restored[column])),column
        assert np.array_equal(restored['master_id'],after['master_id'][selected])
        assert np.all(restored['n_added_bands']==[len(x.split(',')) if x else 0 for x in restored['added_bands']])
        if kind=='NEW_FULL_AVAILABLE':
            assert np.all(result['n_7dt_bands']==n_active)
            if field==1: assert len(result)==466
        else:
            assert np.all(result['is_entirely_new_7dt_source'])
            if field==1: assert len(result)==231
            assert not np.any(result['baseline_finalcat_had_7dt'])
            assert np.all(np.asarray(result['master_class']).astype(str)=='matched')
        print(path,'rows',len(result),'entirely new 7DT',int(missing_parent.sum()),
              'multiple old-parent candidates',int(np.sum(result['baseline_n_source_candidates']>1)),flush=True)
        outputs.append(dict(COSMOS=field,available_bands=n_active,population=kind,N_sources=len(result),
                            N_entirely_new_7dt=int(missing_parent.sum()),filename=name))
    return outputs


def main(fields=(1,2,3,4)):
    """Export the two subsets for each requested field; do not add directories."""
    rows=[]
    for field in fields:
        rows.extend(export_field(field))
    summary=Table(rows=rows)
    print(summary)
    return summary


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fields',type=int,nargs='+',choices=(1,2,3,4),default=[1,2,3,4])
    args=parser.parse_args()
    main(args.fields)
