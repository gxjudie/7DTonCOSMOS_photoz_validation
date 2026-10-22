"""Export three COSMOS1 examples using the existing plot_source_info method.

These illustrative clean SDSS-galaxy examples span 19->20, 17->20 and 16->20
bands and aperture flags 2, 16 and 18. They are not a random statistical sample.
Match reductions by the unique FINALCAT row index, then require every retained
band to use the identical raw detection. Local master IDs change across runs.
Only existing images/catalogs are read; no extraction or fitting is performed.
"""
from pathlib import Path
import importlib.util
import json
import warnings
import numpy as np
from astropy.table import Table
from astropy.units import UnitsWarning
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT/'output/experiments/aperture_flag_ablation'
OUT = ROOT/'output/pdf/cosmos1_new_full20_sources'
SOURCE_IDS = (350, 317, 764)  # Current COSMOS1 field-local master IDs.
BANDS = tuple(f'm{w}' for w in range(400,900,25))


def main():
    warnings.filterwarnings('ignore',category=UnitsWarning)
    before = Table.read(EXP/'baseline/output/COSMOS1/singlephot/cosmos1_master_catalog.fits')
    after = Table.read(ROOT/'output/COSMOS1/singlephot/cosmos1_master_catalog.fits')
    flags_table = Table.read(EXP/'statistical_analysis/cosmos1_new_strict_full20_galaxy_flags.csv')
    spec = importlib.util.spec_from_file_location('sourceinfo_report',ROOT/'Notebooks/sourceinfo_report.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    OUT.mkdir(parents=True,exist_ok=True)
    manifest, band_rows = [], []
    for mid in SOURCE_IDS:
        current = after[np.asarray(after['master_id'])==mid]
        assert len(current)==1
        row = current[0]
        reference_idx = int(row['finalcat_idx'])
        previous = before[(np.asarray(np.ma.filled(before['finalcat_idx'],-1))==reference_idx)
                          & np.asarray(before['has_7dt'],bool)]
        assert len(previous)==1
        old = previous[0]
        flag_rows = flags_table[np.asarray(flags_table['branch_master_id'])==mid]
        assert len(flag_rows)==1
        flags = flag_rows[0]
        added = []
        for band in BANDS:
            sid = row[band+'_source_id']
            assert not np.ma.is_masked(sid) and int(sid)>=0
            assert np.isfinite(row[band+'_mag']) and row[band+'_mag_err']>0
            previous_sid = old[band+'_source_id']
            was_present = not np.ma.is_masked(previous_sid) and int(previous_sid)>=0
            if was_present:
                assert int(previous_sid)==int(sid), 'Retained-band association changed'
            else:
                added.append(band)
                flag = int(flags[band+'_flag'])
                assert flag!=0
                band_rows.append(dict(current_cosmos1_master_id=mid,
                    baseline_cosmos1_master_id=int(old['master_id']),finalcat_idx=reference_idx,
                    band=band,aperture_flag=flag,source_id=int(sid),
                    current_mag=float(row[band+'_mag']),current_mag_err=float(row[band+'_mag_err'])))
        assert int(old['n_7dt_bands'])==20-len(added)<20
        result = module.plot_source_info(after,source_id=mid,field=1,
            project_root=ROOT,pixel_range=40,show=False,
            save_pdf=OUT/f'cosmos1_master{mid}_new_full20_sourceinfo.pdf')
        assert len(result['report_figures'])==2
        item=dict(current_cosmos1_master_id=mid,baseline_cosmos1_master_id=int(old['master_id']),
                  finalcat_idx=reference_idx,sdss_objid=int(row['sdss_objid']),
                  RA_deg=float(row['master_ra']),DEC_deg=float(row['master_dec']),
                  baseline_n_bands=int(old['n_7dt_bands']),current_n_bands=20,
                  added_bands=added,added_band_flags={b:int(flags[b+'_flag']) for b in added},
                  pdf_path=str(result['pdf_path']))
        manifest.append(item)
        print(json.dumps(item),flush=True)
    Table(rows=band_rows).write(OUT/'new_band_details.csv',format='ascii.csv',overwrite=True)
    (OUT/'selection_manifest.json').write_text(json.dumps(dict(
        selection='Illustrative clean one-to-one SDSS galaxies, spanning 1/3/4 recovered bands; not random.',
        identifiers='COSMOS1 field-local master IDs, not combined-catalog global IDs.',
        retained_raw_associations_unchanged=True,photometry_rerun=False,sources=manifest),indent=2)+'\n')


if __name__=='__main__':
    main()
