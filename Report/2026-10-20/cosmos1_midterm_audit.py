"""Compare the June midterm sample with current COSMOS1 photometry.

Run with the existing photoz environment; no production products are changed:
    /opt/anaconda3/envs/photoz/bin/python Report/2026-10-20/cosmos1_midterm_audit.py

This is a photometry and selection audit, not a new photo-z measurement.
Old and new integer source IDs are unrelated. Join catalogs by sky position,
retain the old IDs only as report labels, and record every join separation.
"""
from pathlib import Path
import json
import warnings

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings('ignore', category=u.UnitsWarning)
ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = Path('/Users/mac_seo/Library/Mobile Documents/com~apple~CloudDocs/Study/Research/Undergraduate_internship/7DTonCOSMOS')
OLD_DIR = ARCHIVE / 'Preious_results/2602_results_of_combined_methods'
OUT = ROOT / 'output/COSMOS1/singlephot/analysis/midterm_audit_20261020'
BANDS = tuple(f'm{x}' for x in range(400, 900, 25))
# Values transcribed from Table 4 of the supplied report, not refitted here.
REPORT_IDS = np.array([1853, 1901, 1908, 2122, 2595, 2729, 3833, 3968, 4118, 4236])
REPORT_SPEC_Z = np.array([.166076, .217795, .069669, .125592, .126735, .124199, .185325, .092608, .125197, .099492])
REPORT_PHOTO_Z = np.array([.165, .210, .080, .130, .125, .110, .170, .085, .140, .100])


def numbers(table, name):
    """Convert masked FITS columns to NaN, never to zero or a sentinel flux."""
    return np.asarray(np.ma.filled(table[name], np.nan), dtype=float)


def coords(table, ra, dec):
    return SkyCoord(numbers(table, ra) * u.deg, numbers(table, dec) * u.deg)


def write_csv(rows, filename):
    Table(rows=rows).write(OUT / filename, format='ascii.csv', overwrite=True)


def cigale_input(mags, errors, filename):
    """Convert AB magnitudes directly to mJy; absent measurements remain NaN.

    The old catalog's *_flux columns are f_lambda, while current raw *_flux
    columns are instrumental aperture counts. Neither is directly CIGALE mJy.
    Direct conversion from calibrated AB magnitude avoids that unit ambiguity.
    Negative input redshift enables photo-z mode, but a model/redshift grid
    must still be supplied before any CIGALE fit is meaningful.
    """
    t = Table({'id': REPORT_IDS, 'redshift': np.full(10, -1)})
    for bi, band in enumerate(BANDS):
        flux = 10 ** ((23.9 - mags[:, bi]) / 2.5) / 1000
        t[f'7DT_{band}'] = flux
        t[f'7DT_{band}_err'] = flux * np.log(10) / 2.5 * errors[:, bi]
    t.write(OUT / filename, format='ascii.basic', overwrite=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old_path = OLD_DIR / 'COSMOS1_master_catalog_with_SDSS.csv'
    long_path = OLD_DIR / 'COSMOS1_cat_combined_has20_complete.csv'
    current_path = ROOT / 'output/COSMOS1/singlephot/cosmos1_master_catalog.fits'
    old = Table.read(old_path, format='ascii.csv')
    old_long = Table.read(long_path, format='ascii.csv')
    current = Table.read(current_path)
    current = current[np.asarray(current['has_7dt'], dtype=bool)]
    oc = coords(old, 'ra', 'dec')
    nc = coords(current, 'ra_7dt', 'dec_7dt')
    nearest, sep, _ = oc.match_to_catalog_sky(nc)
    _, second_sep, _ = oc.match_to_catalog_sky(nc, nthneighbor=2)
    joined = sep.arcsec < .7
    # Require a unique target for every accepted master-catalog association.
    assert len(np.unique(nearest[joined])) == joined.sum()
    old_index = np.array([np.flatnonzero(old['id'] == oid).item() for oid in REPORT_IDS])
    new_index = nearest[old_index]
    assert np.all(joined[old_index])
    assert np.all(second_sep.arcsec[old_index] > .7)
    old_mag = np.column_stack([numbers(old, f'{b}_mag')[old_index] for b in BANDS])
    old_err = np.column_stack([numbers(old, f'{b}_mag_err')[old_index] for b in BANDS])
    new_mag = np.column_stack([numbers(current, f'{b}_mag')[new_index] for b in BANDS])
    new_err = np.column_stack([numbers(current, f'{b}_mag_err')[new_index] for b in BANDS])
    raw_mag = np.full_like(old_mag, np.nan)
    raw_err = np.full_like(old_err, np.nan)
    raw_flags = np.zeros_like(old_mag, dtype=int)
    rows, band_stats = [], []
    for bi, band in enumerate(BANDS):
        raw_path = ROOT / f'output/COSMOS1/singlephot/cat_cosmos1_{band}.fits'
        raw = Table.read(raw_path)
        rc = coords(raw, 'world_centroid_alpha', 'world_centroid_delta')
        # Compare each old band centroid to its new raw-band neighbor. This
        # differs from using the average master position for every filter.
        selected_old = old_long[np.asarray(old_long['band']) == band]
        lookup = {int(r['id']): r for r in selected_old}
        selected_old = Table(rows=[lookup[int(oid)] for oid in REPORT_IDS], names=selected_old.colnames)
        bc = coords(selected_old, 'world_centroid_alpha', 'world_centroid_delta')
        ri, rsep, _ = bc.match_to_catalog_sky(rc)
        _, second_raw, _ = bc.match_to_catalog_sky(rc, nthneighbor=2)
        # A generous diagnostic radius is explicit. It does not override the
        # production quality selection, nor establish that flagged flux is safe.
        assert np.all(rsep.arcsec < 1.5)
        raw_mag[:, bi] = numbers(raw, 'aper3_mag')[ri]
        raw_err[:, bi] = numbers(raw, 'aper3_mag_err')[ri]
        raw_flags[:, bi] = np.asarray(raw['aper3_flags'][ri], dtype=int)
        for k, oid in enumerate(REPORT_IDS):
            present = bool(np.isfinite(new_mag[k, bi]))
            if present:
                assert current[f'{band}_source_id'][new_index[k]] == raw['source_id'][ri[k]]
            else:
                assert raw_flags[k, bi] != 0
            rows.append(dict(old_id=int(oid), current_master_id=int(current['master_id'][new_index[k]]),
                             band=band, old_mag=float(old_mag[k, bi]), current_master_mag=float(new_mag[k, bi]),
                             current_raw_mag=float(raw_mag[k, bi]), old_mag_err=float(old_err[k, bi]),
                             current_master_mag_err=float(new_err[k, bi]), current_raw_mag_err=float(raw_err[k, bi]),
                             old_aperture_flags=int(selected_old['aperture_flags'][k]),
                             current_aperture_flags=int(raw_flags[k, bi]), included_in_current_master=present,
                             raw_source_id=int(raw['source_id'][ri[k]]), old_to_raw_sep_arcsec=float(rsep.arcsec[k]),
                             second_raw_sep_arcsec=float(second_raw.arcsec[k]),
                             raw_match_crowded=bool(second_raw.arcsec[k] < 2*rsep.arcsec[k]),
                             delta_raw_mag=float(raw_mag[k, bi]-old_mag[k, bi])))
        dm = numbers(current, f'{band}_mag')[nearest] - numbers(old, f'{band}_mag')
        usable = joined & np.isfinite(dm)
        band_stats.append(dict(band=band, n_matched_measurements=int(usable.sum()),
                               median_current_minus_old_mag=float(np.median(dm[usable])),
                               p16_delta_mag=float(np.percentile(dm[usable], 16)),
                               p84_delta_mag=float(np.percentile(dm[usable], 84))))
    delta = new_mag-old_mag
    raw_delta = raw_mag-old_mag
    included = np.isfinite(new_mag)
    summary_rows = []
    for k, oid in enumerate(REPORT_IDS):
        summary_rows.append(dict(old_id=int(oid), current_master_id=int(current['master_id'][new_index[k]]),
                                 master_sep_arcsec=float(sep.arcsec[old_index[k]]),
                                 z_spec_report=float(REPORT_SPEC_Z[k]), z_phot_report=float(REPORT_PHOTO_Z[k]),
                                 old_n_bands=20, current_n_bands=int(included[k].sum()),
                                 dropped_nonzero_flag_bands=','.join(np.array(BANDS)[~included[k]]),
                                 median_delta_mag_retained=float(np.nanmedian(delta[k])),
                                 std_delta_mag_retained=float(np.nanstd(delta[k])),
                                 median_error_ratio_retained=float(np.nanmedian(new_err[k]/old_err[k])),
                                 max_abs_delta_mag_raw=float(np.max(np.abs(raw_delta[k])))))
    write_csv(rows, 'ten_galaxies_band_comparison.csv')
    write_csv(summary_rows, 'ten_galaxies_summary.csv')
    write_csv(band_stats, 'old_parent_per_band_comparison.csv')
    cigale_input(old_mag, old_err, 'cigale_old_ten.txt')
    cigale_input(new_mag, new_err, 'cigale_current_master_ten.txt')
    cigale_input(raw_mag, raw_err, 'cigale_current_raw_ten_DIAGNOSTIC_ONLY.txt')

    all_mag = np.column_stack([numbers(current, f'{b}_mag') for b in BANDS])
    all_err = np.column_stack([numbers(current, f'{b}_mag_err') for b in BANDS])
    success = np.isfinite(all_mag) & (all_mag > 0) & (all_mag < 50) & np.isfinite(all_err) & (all_err > 0)
    nband = success.sum(1)
    z = numbers(current, 'finalcat_Z_TOT_Z')
    source = np.char.strip(np.asarray(current['finalcat_Z_TOT_ZSOURCE']).astype(str))
    # This is a redshift-availability cut, not a spectroscopic quality cut.
    spec = np.isfinite(z) & (z > 0) & (source != 'NN')
    groups = [('current_7dt', np.ones(len(current), bool)), ('current_7dt_with_reference_z', spec)]
    pop_rows = []
    for label, use in groups:
        for nb in range(21):
            pop_rows.append(dict(population=label, n_bands=nb, n_sources=int(np.sum(use & (nband == nb)))))
    write_csv(pop_rows, 'current_population_band_counts.csv')
    dz = REPORT_SPEC_Z-REPORT_PHOTO_Z
    stats = dict(old_master_path=str(old_path), old_long_path=str(long_path), current_master_path=str(current_path),
                 old_n_sources=len(old), old_matched_within_0p7_arcsec=int(joined.sum()),
                 old_parent_matched_and_current_20bands=int(np.sum(joined & (nband[nearest] == 20))),
                 ten_current_band_counts=included.sum(1).tolist(), ten_master_max_sep_arcsec=float(sep.arcsec[old_index].max()),
                 ten_raw_max_sep_arcsec=max(r['old_to_raw_sep_arcsec'] for r in rows),
                 ten_raw_min_second_neighbor_arcsec=min(r['second_raw_sep_arcsec'] for r in rows),
                 ten_raw_crowded_associations=sum(r['raw_match_crowded'] for r in rows),
                 n_old_measurements=200, n_current_retained=int(included.sum()), n_dropped=int((~included).sum()),
                 dropped_flags={str(flag): int(np.sum(raw_flags[~included] == flag)) for flag in np.unique(raw_flags[~included])},
                 retained_abs_delta_mag_quantiles=np.percentile(np.abs(delta[included]), [50, 90, 95, 100]).tolist(),
                 raw_abs_delta_mag_quantiles=np.percentile(np.abs(raw_delta), [50, 90, 95, 100]).tolist(),
                 retained_median_error_ratio=float(np.nanmedian(new_err/old_err)),
                 report_table_recomputed_metrics=dict(eta=float(np.mean(np.abs(dz)/(1+REPORT_SPEC_Z)>.15)),
                 sigma_nmad=float(1.48*np.median(np.abs(dz-np.median(dz))/(1+REPORT_SPEC_Z))),
                 bias=float(np.mean(-dz/(1+REPORT_SPEC_Z)))),
                 new_photoz_fitted=False,
                 weight_map_test='Attempted isolated m400 rerun with weight maps; installed SourceXtractor++ 1.1.0 crashed before producing a catalog in sepp and sepp_clean environments.',
                 limitation='Senior code/results and original executed CIGALE grid are unavailable. Raw flagged photometry is diagnostic only; five nearest-raw associations have a second candidate within twice the nearest distance.')
    (OUT/'summary.json').write_text(json.dumps(stats, indent=2)+'\n')

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), layout='constrained')
    wave = np.array([int(b[1:]) for b in BANDS])
    for k, oid in enumerate(REPORT_IDS):
        axes[0].plot(wave, delta[k], '.-', alpha=.7, lw=.8, label=str(oid))
    axes[0].axhline(0, color='black', lw=.7)
    axes[0].set(xlabel='Filter center [nm]', ylabel='Current - old magnitude', title='Same ten galaxies: retained bands')
    axes[0].legend(ncol=2, fontsize=7)
    axes[0].grid(alpha=.2)
    im = axes[1].imshow(raw_flags != 0, cmap='RdYlGn_r', vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    axes[1].set(xticks=np.arange(0,20,2), xticklabels=np.array(BANDS)[::2], yticks=np.arange(10),
                yticklabels=[f'{oid}: {n}/20' for oid,n in zip(REPORT_IDS,included.sum(1))], title='Red: raw flag != 0, dropped from master')
    axes[1].tick_params(axis='x', rotation=60)
    for k in range(10):
        for bi in range(20):
            if raw_flags[k,bi]: axes[1].text(bi,k,str(raw_flags[k,bi]),ha='center',va='center',fontsize=6)
    for label, use in [('All 7DT', np.ones(len(current),bool)), ('Reference z available',spec)]:
        counts=np.bincount(nband[use],minlength=21)
        axes[2].step(np.arange(21), counts, where='mid',label=f'{label} (N={use.sum():,})')
    axes[2].set(xlabel='Successful bands in current COSMOS1 master',ylabel='Sources',yscale='log',
                xticks=[0,5,10,15,20],title='Current parent population differs')
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=.2)
    fig.suptitle('COSMOS1 midterm audit: photometry and selection; no new photo-z fit',fontsize=13)
    fig.savefig(OUT/'cosmos1_midterm_audit.png',dpi=180)
    plt.close(fig)
    # Show the SED changes directly. Flagged raw measurements are crosses,
    # not reinstated science measurements, especially in crowded detections.
    fig, axes = plt.subplots(2, 5, figsize=(17, 7), layout='constrained')
    for k, (ax, oid) in enumerate(zip(axes.ravel(), REPORT_IDS)):
        ax.errorbar(wave, old_mag[k], yerr=old_err[k], fmt='o-', ms=3,
                    lw=.8, color='0.5', label='Old: 20 bands')
        ax.errorbar(wave, new_mag[k], yerr=new_err[k], fmt='o', ms=3,
                    color='tab:blue', label='Current master')
        dropped = ~included[k]
        ax.scatter(wave[dropped], raw_mag[k,dropped], marker='x',
                   color='tab:red', s=32, label='Flagged raw candidate')
        ax.invert_yaxis()
        ax.set(title=f'Old ID {oid}: now {included[k].sum()}/20',
               xlabel='Filter center [nm]', ylabel='AB aperture magnitude')
        ax.grid(alpha=.2)
    axes.ravel()[0].legend(fontsize=7)
    fig.suptitle('Same ten midterm galaxies: old and current COSMOS1 SEDs', fontsize=14)
    fig.savefig(OUT/'ten_galaxies_sed_comparison.png',dpi=160)
    plt.close(fig)
    print(json.dumps(stats,indent=2))
    print('OUTPUT',OUT)


if __name__ == '__main__':
    main()
