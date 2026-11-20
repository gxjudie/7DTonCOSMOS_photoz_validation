"""Prepare compact baseline presentation assets and audit sky/source offsets.

Read existing blank-aperture grids and the agreed comparison CSV. Never change
the production catalogs or the CSV's science definitions. The 20-band overview
uses actual COSMOS IDs and the deepest valid field at each overlapping position.
"""

from pathlib import Path
import csv
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

warnings.filterwarnings('ignore', category=FITSFixedWarning)
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
BASE = ROOT / 'output/experiments/aperture_flag_ablation/baseline/full_products'
BANDS = tuple(f'm{w}' for w in range(400, 900, 25))
FIELDS = (1, 2, 3, 4)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12})


def numbers(column):
    """Restore missing FITS values as NaN, including integer columns."""
    return np.asarray(np.ma.filled(column.astype(float), np.nan), float)


def labels(column):
    return np.char.strip(np.asarray(column).astype(str))


sky_tables, masters, wcs_by_field, shapes = {}, {}, {}, {}
all_corners = []
for field in FIELDS:
    product = BASE / f'COSMOS{field}/singlephot'
    sky_tables[field] = Table.read(product / 'analysis' / f'cosmos{field}_7dt_depth_5sigma_3arcsec.csv')
    masters[field] = Table.read(product / f'cosmos{field}_master_catalog.fits')
    path = next((ROOT / f'COSMOS_data/COSMOS{field}').glob('*m400*.com.fits'))
    header = fits.getheader(path)
    wcs_by_field[field] = WCS(header)
    nx, ny = header['NAXIS1'], header['NAXIS2']
    shapes[field] = (ny, nx)
    # Pixel edges are in FITS-origin coordinates: zero-based array centers are
    # 1..N in WCS origin=1. Match the archived map's grid-edge convention below.
    all_corners.append(wcs_by_field[field].all_pix2world([[0, 0], [nx, 0], [nx, ny], [0, ny]], 1))

corners = np.concatenate(all_corners)
ra_low, ra_high = np.min(corners[:, 0]), np.max(corners[:, 0])
dec_low, dec_high = np.min(corners[:, 1]), np.max(corners[:, 1])
world_ra, world_dec = np.meshgrid(np.linspace(ra_low, ra_high, 500), np.linspace(dec_low, dec_high, 350))
lookup = {}
for field in FIELDS:
    ny, nx = shapes[field]
    x, y = wcs_by_field[field].all_world2pix(world_ra, world_dec, 1)
    ix = np.searchsorted(np.linspace(0, nx, 11), x, side='right') - 1
    iy = np.searchsorted(np.linspace(0, ny, 11), y, side='right') - 1
    inside = (ix >= 0) & (ix < 10) & (iy >= 0) & (iy < 10)
    lookup[field] = (ix, iy, inside)


def grid_for(field, band):
    """Keep only originally valid sky cells; unavailable bands remain all NaN."""
    table = sky_tables[field]
    grid = np.full((10, 10), np.nan)
    selected = (labels(table['band']) == band) & (labels(table['status']) == 'ok')
    for row in table[selected]:
        value = float(row['depth_5sigma_mag'])
        if np.isfinite(value):
            grid[int(row['grid_row']), int(row['grid_col'])] = value
    return grid


# Redraw scientific depth data with one common colorbar, rather than shrinking
# twenty existing screenshots with repeated axis labels and colorbars.
fig, axes = plt.subplots(4, 5, figsize=(18, 10.4), sharex=True, sharey=True)
fig.subplots_adjust(left=0.055, right=0.91, bottom=0.11, top=0.92, wspace=0.08, hspace=0.23)
cmap = plt.get_cmap('viridis').copy()
cmap.set_bad('#E5E5E5')
for band, ax in zip(BANDS, axes.flat):
    combined = np.full(world_ra.shape, np.nan)
    for field in FIELDS:
        grid = grid_for(field, band)
        ix, iy, inside = lookup[field]
        candidate = np.full(world_ra.shape, np.nan)
        candidate[inside] = grid[iy[inside], ix[inside]]
        combined = np.fmax(combined, candidate)
    image = ax.imshow(np.ma.masked_invalid(combined), origin='lower',
                      extent=[ra_low, ra_high, dec_low, dec_high],
                      cmap=cmap, norm=Normalize(13, 22.5), interpolation='nearest', aspect='auto')
    # Explicit limits avoid toggling a shared x axis once per subplot.
    ax.set_xlim(ra_high, ra_low)
    ax.set_title(band, fontsize=15, pad=4)
    ax.set_xticks([149.5, 150.0, 150.5, 151.0])
    ax.set_yticks([1.5, 2.0, 2.5])
    ax.tick_params(labelsize=9, length=3)
    ax.spines[['top', 'right']].set_visible(False)
for ax in axes[-1]:
    ax.set_xlabel('RA [deg]', fontsize=11)
for ax in axes[:, 0]:
    ax.set_ylabel('Dec [deg]', fontsize=11)
color_ax = fig.add_axes([0.93, 0.18, 0.014, 0.66])
fig.colorbar(image, cax=color_ax, label='5σ limiting magnitude [AB mag]')
fig.suptitle('7DT COSMOS1–4: 20-band sky depth maps', fontsize=22, y=0.98)
fig.text(0.055, 0.045, '3 arcsec aperture. Baseline: aperture flag = 0. Overlaps show the deepest valid field.', fontsize=12)
fig.text(0.055, 0.018, 'Gray: no valid depth estimate. All fields shown, including the COSMOS3-only region.', fontsize=11, color='#444444')
fig.savefig(OUT / 'cosmos1-4_20band_depth_maps_baseline.png', dpi=220, facecolor='white')
plt.close(fig)


def snr_locus(field, band):
    """Reproduce the baseline matched-sample S/N and 0.25-mag bin medians."""
    master = masters[field]
    mag, err = numbers(master[f'{band}_mag']), numbers(master[f'{band}_mag_err'])
    with np.errstate(divide='ignore', invalid='ignore'):
        snr = 1.0857 / err
    selected = (labels(master['master_class']) == 'matched') & np.isfinite(mag) & np.isfinite(snr) & (err > 0)
    values, sn = mag[selected], snr[selected]
    edges = np.arange(np.floor(values.min()), np.ceil(values.max()) + 0.25, 0.25)
    centers, medians = [], []
    for left, right in zip(edges[:-1], edges[1:]):
        use = (values >= left) & (values < right)
        if use.sum() >= 10:
            centers.append((left + right) / 2)
            medians.append(float(np.median(sn[use])))
    return mag, snr, selected, np.asarray(centers), np.asarray(medians)


with (OUT / 'depth_comparison_baseline.csv').open() as stream:
    comparison = {r['band']: r for r in csv.DictReader(stream)}

# Show both an ordinary example and the largest-offset comparison, avoiding
# selection of only the closest-agreement bands in the presentation narrative.
fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.3))
for ax, (field, band) in zip(axes, [(1, 'm600'), (3, 'm475')]):
    mag, snr, selected, centers, medians = snr_locus(field, band)
    source = float(comparison[band][f'source5snr_cosmos{field}_mag'])
    ax.scatter(mag[selected], snr[selected], s=8, alpha=0.25, color='#2878B5', rasterized=True)
    ax.plot(centers, medians, 'o-', color='#C63932', markersize=4, linewidth=1.5, label='Bin median (N ≥ 10)')
    ax.axhline(5, color='black', linestyle='--', linewidth=1.2, label='S/N = 5')
    ax.axvline(source, color='#D18A21', linestyle='-.', linewidth=1.4, label=f'5σ: {source:.2f} mag')
    ax.set(xlim=(source - 1.0, source + 0.65), ylim=(0, 10), xlabel=f'7DT {band} magnitude [AB]', ylabel='S/N')
    ax.set_title(f'COSMOS{field} {band}', fontsize=17)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=10)
fig.text(0.07, 0.02, 'Baseline matched sources. S/N = 1.0857 / magnitude error. 0.25-mag bins with at least 10 sources.', fontsize=10)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(OUT / 'source_depth_examples_baseline.png', dpi=200, facecolor='white')
plt.close(fig)


# Assess ALL estimable field-band pairs, not only the illustrative bars. Sky
# sigma is spatial scatter and must never be treated as an estimator error bar
# or used to assign a statistical significance to the method difference.
audit = []
for band in BANDS:
    for field in FIELDS:
        row = comparison[band]
        sky = float(row[f'sky5snr_cosmos{field}_mag_avg'])
        sigma = float(row[f'sky5snr_cosmos{field}_mag_sigma'])
        source = float(row[f'source5snr_cosmos{field}_mag'])
        if not (np.isfinite(sky) and np.isfinite(source)):
            continue
        master = masters[field]
        mag, snr, selected, _, _ = snr_locus(field, band)
        x, y = numbers(master['master_x']) - 1, numbers(master['master_y']) - 1
        ny, nx = shapes[field]
        near = selected & (snr >= 4) & (snr <= 6) & np.isfinite(x) & np.isfinite(y)
        indices = np.flatnonzero(near)
        grid = grid_for(field, band)
        local_sky, source_noise_limits = [], []
        for index in indices:
            gx, gy = int(np.floor(x[index] / nx * 10)), int(np.floor(y[index] / ny * 10))
            if 0 <= gx < 10 and 0 <= gy < 10 and np.isfinite(grid[gy, gx]):
                local_sky.append(grid[gy, gx])
                # Translate each near-threshold source's flux/error to its own
                # 5-sigma magnitude. This is a diagnostic noise-limit proxy,
                # not a replacement for the agreed median-crossing estimator.
                source_noise_limits.append(mag[index] + 2.5 * np.log10(snr[index] / 5))
        local_sky, source_noise_limits = np.asarray(local_sky), np.asarray(source_noise_limits)
        delta = source - sky
        audit.append(dict(field=field, band=band, sky_mean_mag=sky, sky_spatial_sigma_mag=sigma,
                          source_mag=source, delta_source_minus_sky_mag=delta,
                          flux_limit_ratio_sky_over_source=10 ** (0.4 * delta),
                          n_local_near_snr5=len(local_sky),
                          median_local_sky_mag=float(np.median(local_sky)) if len(local_sky) else np.nan,
                          median_source_noise_limit_mag=float(np.median(source_noise_limits)) if len(local_sky) else np.nan,
                          median_local_offset_mag=float(np.median(source_noise_limits-local_sky)) if len(local_sky) else np.nan))
with (OUT / 'sky_source_offset_audit.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(audit[0]))
    writer.writeheader()
    writer.writerows(audit)

fig, ax = plt.subplots(figsize=(11.8, 4.5))
field_colors = ['#2878B5', '#D18A21', '#C63932', '#348757']
for field, color in zip(FIELDS, field_colors):
    rows = [r for r in audit if r['field'] == field]
    ax.scatter([int(r['band'][1:]) for r in rows], [r['delta_source_minus_sky_mag'] for r in rows],
               label=f'COSMOS{field}', color=color, s=36, zorder=3)
ax.axhline(0, color='black', linestyle='--', linewidth=1)
ax.set(xlabel='Band central wavelength [nm]', ylabel='Source depth − sky depth [mag]',
       ylim=(-0.03, 0.53), xlim=(385, 890))
ax.set_title('Sky/source offsets across all estimable field-band pairs', fontsize=16)
ax.legend(ncol=4, frameon=False, fontsize=11)
ax.grid(alpha=0.2)
fig.text(0.09, 0.02, 'Positive: source-based depth is deeper. Missing estimates are omitted, not treated as agreement.', fontsize=10)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(OUT / 'sky_source_offset_audit.png', dpi=200, facecolor='white')
plt.close(fig)

delta = np.array([r['delta_source_minus_sky_mag'] for r in audit])
print('Compared pairs:', len(audit))
print('Offset min/median/max:', delta.min(), np.median(delta), delta.max())
print('Local near-SNR5 offsets:')
for row in sorted(audit, key=lambda r: r['delta_source_minus_sky_mag'], reverse=True)[:5]:
    print(row)
print('Saved the 20-band overview, source examples, and separate offset audit.')
