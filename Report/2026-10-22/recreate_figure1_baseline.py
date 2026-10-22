"""Recreate proposal Figure 1 from the saved pre-aperture-flag-ablation results.

Run with /opt/anaconda3/envs/photoz/bin/python. This standalone copy of the
archived Proposal.ipynb plotting cells deliberately reads the baseline catalog
and baseline depth CSVs, never the current no-flag-cut production products.
No photometry, matching, catalog selection, or branch state is changed.
The only display changes are a white/black grayscale scheme for both depth panels,
600-dpi PNG output, and embedded TrueType fonts in the vector PDF.
"""
import hashlib
import json
import matplotlib
matplotlib.use('Agg')

# Independently runnable: read the persisted catalog, without rebuilding it
# or repeating the SDSS queries in the catalog-production section.
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from astropy.io import fits
from astropy.table import Table
from scl_plot import scl_plt

CUM_ROOT = Path(__file__).resolve().parents[2]
BASELINE = CUM_ROOT / 'output/experiments/aperture_flag_ablation/baseline/full_products'
FINAL_OUTPUT = CUM_ROOT / 'Report/7DT Observation Proposal/Figure1_2026-10-22_baseline_flagcut'
FINAL_OUTPUT.mkdir(parents=True, exist_ok=True)
CUM_OUTPUT = FINAL_OUTPUT
# The archive retains aper3_flags == 0 and the original success masks.
CUM_PATH = BASELINE / 'singlephot/master_catalog_MatchedWithFianlcat_OverallCosmosField.fits'
CUM_CATALOG = Table.read(CUM_PATH, hdu='CATALOG')
CUM_REFERENCE = Table.read(CUM_ROOT / 'COSMOS_data/COSMOS_finalcat.fits')
CUM_BANDS = tuple(f'm{wave}' for wave in range(400, 900, 25))
CUM_FIELDS = (1, 2, 3, 4)
with fits.open(CUM_PATH, memmap=False) as hdus:
    # The production file stores exact masks separately, including Boolean
    # columns. Honor those masks when reading the saved per-band success flags.
    CUM_PACKED_MASKS = np.array(hdus['COLUMN_MASKS'].data, copy=True)
    CUM_BITORDER = hdus['COLUMN_MASKS'].header['BITORDER']
CUM_COLUMN_INDEX = {name: i for i, name in enumerate(CUM_CATALOG.colnames)}


def cumulative_catalog_values(name, fill=np.nan):
    """Return one saved column with its production FITS mask restored."""
    index = CUM_COLUMN_INDEX[name]
    bit = index % 8 if CUM_BITORDER == 'little' else 7 - index % 8
    masked = ((CUM_PACKED_MASKS[:, index // 8] >> bit) & 1).astype(bool)
    values = np.asarray(CUM_CATALOG[name], dtype=float).copy()
    values[masked | np.ma.getmaskarray(CUM_CATALOG[name])] = fill
    return values


CUM_RM = cumulative_catalog_values('finalcat_RM')
CUM_REF_RM = np.asarray(np.ma.filled(CUM_REFERENCE['RM'], np.nan), dtype=float)
CUM_VALID = np.isfinite(CUM_RM) & (CUM_RM > 0) & (CUM_RM < 50)
CUM_REF_VALID = np.isfinite(CUM_REF_RM) & (CUM_REF_RM > 0) & (CUM_REF_RM < 50)
# Keep the previous section's 0.5-mag spacing over the full valid RM range.
CUM_ALL_RM = np.concatenate([CUM_RM[CUM_VALID], CUM_REF_RM[CUM_REF_VALID]])
CUM_THRESHOLDS = np.arange(np.floor(CUM_ALL_RM.min() * 2) / 2,
                           np.ceil(CUM_ALL_RM.max() * 2) / 2 + 0.5, 0.5)


def cumulative_rm_counts(values):
    """Count RM <= each threshold, including values on the threshold itself.

    Count from the complete valid bright-end sample before restricting the
    display to 16-21 mag; a cropped histogram would omit brighter sources.
    """
    return np.searchsorted(np.sort(values), CUM_THRESHOLDS, side='right')


CUM_REFERENCE_COUNTS = cumulative_rm_counts(CUM_REF_RM[CUM_REF_VALID])
CUM_BAND_COUNTS = {}
for band in CUM_BANDS:
    successful = np.zeros(len(CUM_CATALOG), dtype=bool)
    for field in CUM_FIELDS:
        prefix = f'7dt_cosmos{field}_'
        matched = cumulative_catalog_values(prefix + 'matched', fill=0).astype(bool)
        success = cumulative_catalog_values(prefix + band + '_success', fill=0).astype(bool)
        # The saved success flags encode source_id >= 0, 0 < mag < 50,
        # finite mag/error, and mag_err > 0. Union fields within each band.
        successful |= matched & success
    CUM_BAND_COUNTS[band] = cumulative_rm_counts(CUM_RM[successful & CUM_VALID])

CUM_TABLE = Table([CUM_THRESHOLDS, CUM_REFERENCE_COUNTS] +
                  [CUM_BAND_COUNTS[band] for band in CUM_BANDS],
                  names=['RM_threshold', 'N_reference'] + [f'N_{band}' for band in CUM_BANDS])
CUM_WAVELENGTH = np.array([int(band[1:]) for band in CUM_BANDS], dtype=float)
CUM_CMAP = LinearSegmentedColormap.from_list(
    'cumulative_purple_to_red',
    ['#8A00CC', '#003BFF', '#00BFFF', '#00A83B', '#FFD400', '#FF6500', '#E00000'])
CUM_NORM = Normalize(vmin=CUM_WAVELENGTH.min(), vmax=CUM_WAVELENGTH.max())

# Match the preceding figure: cumulative bin counts at bin centers,
# with connected lines and small circular markers for the 7DT bands.
CUM_BIN_CENTERS = 0.5 * (CUM_THRESHOLDS[:-1] + CUM_THRESHOLDS[1:])



# Compare every threshold and every band with the persisted pre-ablation table.
# This detects accidentally using the promoted no-flag-cut catalog.
BASELINE_COUNTS = Table.read(BASELINE / 'singlephot/analysis/numofmatchedbands/rm_cumulative_source_count_by_band.csv')
assert CUM_TABLE.colnames == BASELINE_COUNTS.colnames
for name in CUM_TABLE.colnames:
    np.testing.assert_array_equal(CUM_TABLE[name], BASELINE_COUNTS[name])
CUM_TABLE.write(FINAL_OUTPUT / 'figure1_baseline_cumulative_counts.csv', format='ascii.csv', overwrite=True)
print('All cumulative source counts equal the archived baseline, in every band.')

# Keep the combined-depth section independently runnable from root/Notebooks.
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

PROJECT_ROOT = CUM_ROOT
DATA_ROOT = PROJECT_ROOT / 'COSMOS_data'
BANDS_7DT = tuple(f'm{wave}' for wave in range(400, 900, 25))
OVERALL_ANALYSIS_DIR = FINAL_OUTPUT


def _combined_depth_band_image_path(cosmos, band):
    """Resolve the band's coadd for the same WCS used by Proposal.ipynb."""
    paths = sorted((DATA_ROOT / f'COSMOS{cosmos}').glob(f'*{band}*.com.fits'))
    if len(paths) != 1:
        raise FileNotFoundError(f'Expected one COSMOS{cosmos} {band} coadd; got {len(paths)}')
    return paths[0]

# Load the COSMOS1-4 depth grids and convert their cell boundaries to RA/Dec.
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon

COMBINED_DEPTH_FIG_DIR = OVERALL_ANALYSIS_DIR
COMBINED_DEPTH_FIG_DIR.mkdir(parents=True, exist_ok=True)
COMBINED_DEPTH_COSMOS = (1, 2, 3, 4)
# Display-only numbering follows mathematical quadrants on the reversed-RA
# plot: upper right=1, upper left=2, lower left=3, lower right=4.
# Original COSMOS IDs still select every input file, WCS, and depth grid.
COMBINED_DEPTH_DISPLAY_FIELD = {1: 4, 2: 1, 3: 3, 4: 2}
COMBINED_DEPTH_POOL = 16
COMBINED_DEPTH_PAD_DEG = 0.02
COMBINED_DEPTH_BACKGROUND_ALPHA = 0.22
COMBINED_DEPTH_LABEL_BBOX_ALPHA = 0.30
COMBINED_DEPTH_MAP_NX = 1000
COMBINED_DEPTH_MAP_NY = 750


def _combined_depth_csv_path(cosmos):
    return (BASELINE / f'COSMOS{cosmos}' / 'singlephot' / 'analysis'
            / f'cosmos{cosmos}_7dt_depth_5sigma_3arcsec.csv')


def _combined_depth_detection_path(cosmos):
    candidates = [
        DATA_ROOT / f'COSMOS{cosmos}' / f'COSMOS{cosmos}_detectionimage.fits',
        PROJECT_ROOT / 'output' / f'COSMOS{cosmos}' / f'COSMOS{cosmos}_detectionimage.fits',
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f'No detection image found for COSMOS{cosmos}')


def _combined_depth_pooled_background(path, factor=COMBINED_DEPTH_POOL):
    data = fits.getdata(path, memmap=True)
    ny, nx = data.shape
    nx_used, ny_used = (nx // factor) * factor, (ny // factor) * factor
    pooled = np.empty((ny_used // factor, nx_used // factor), dtype=np.float32)
    for y0 in range(0, ny_used, 1024):
        y1 = min(y0 + 1024, ny_used)
        chunk = np.asarray(data[y0:y1, :nx_used], dtype=np.float32)
        chunk_h = (y1 - y0) // factor
        block = chunk[:chunk_h * factor].reshape(chunk_h, factor, nx_used // factor, factor)
        with np.errstate(invalid='ignore'):
            pooled[y0 // factor:y0 // factor + chunk_h] = np.nanmean(block, axis=(1, 3))
    return np.ma.masked_equal(pooled, 0.0), nx_used, ny_used


combined_depth_tables = {}
combined_depth_available_by_band = {band: set() for band in BANDS_7DT}
combined_depth_maps = {}
combined_depth_footprints = {}
all_ra, all_dec, background_samples = [], [], []

for cosmos in COMBINED_DEPTH_COSMOS:
    csv_path = _combined_depth_csv_path(cosmos)
    if not csv_path.exists():
        raise FileNotFoundError(f'Missing depth table: {csv_path}')
    table = Table.read(csv_path, format='csv')
    combined_depth_tables[cosmos] = table
    available = set(np.asarray(table['band']).astype(str))
    for band in available:
        combined_depth_available_by_band.setdefault(band, set()).add(cosmos)

    detection_path = _combined_depth_detection_path(cosmos)
    detection_header = fits.getheader(detection_path)
    detection_wcs = WCS(detection_header)
    background, nx_used, ny_used = _combined_depth_pooled_background(detection_path)
    x_corners = np.array([0.5, nx_used + 0.5, nx_used + 0.5, 0.5])
    y_corners = np.array([0.5, 0.5, ny_used + 0.5, ny_used + 0.5])
    ra_corners, dec_corners = detection_wcs.all_pix2world(x_corners, y_corners, 1)
    combined_depth_footprints[cosmos] = {
        'background': background,
        'wcs': detection_wcs,
        'ra_corners': np.asarray(ra_corners, dtype=float),
        'dec_corners': np.asarray(dec_corners, dtype=float),
        'ra_min': float(np.nanmin(ra_corners)),
        'ra_max': float(np.nanmax(ra_corners)),
        'dec_min': float(np.nanmin(dec_corners)),
        'dec_max': float(np.nanmax(dec_corners)),
        'nx_used': nx_used,
        'ny_used': ny_used,
        'polygon': np.column_stack((ra_corners, dec_corners)),
    }
    all_ra.extend(ra_corners.tolist())
    all_dec.extend(dec_corners.tolist())
    finite_background = background.compressed()
    if finite_background.size:
        background_samples.append(finite_background[::max(1, finite_background.size // 50000)])

    for band in available:
        band_path = _combined_depth_band_image_path(cosmos, band)
        band_header = fits.getheader(band_path)
        # Match by the FITS FILTER keyword, never by filename order. Read the
        # coadd's actual EXPTIME: COSMOS1-3 report total exposure duration;
        # COSMOS4 describes its keyword as maximum equivalent exposure time.
        if str(band_header['FILTER']).strip() != band:
            raise ValueError(f'Unexpected FILTER in {band_path}')
        exposure_seconds = float(band_header['EXPTIME'])
        if not np.isfinite(exposure_seconds) or exposure_seconds <= 0:
            raise ValueError(f'Invalid EXPTIME in {band_path}')
        band_wcs = WCS(band_header)
        band_rows = table[np.asarray(table['band']).astype(str) == band]
        x_edges = np.unique(np.concatenate([
            np.asarray(band_rows['x_min'], dtype=float), np.asarray(band_rows['x_max'], dtype=float)
        ]))
        y_edges = np.unique(np.concatenate([
            np.asarray(band_rows['y_min'], dtype=float), np.asarray(band_rows['y_max'], dtype=float)
        ]))
        if len(x_edges) < 2 or len(y_edges) < 2:
            continue
        pixel_x, pixel_y = np.meshgrid(x_edges + 0.5, y_edges + 0.5)
        ra_edges, dec_edges = band_wcs.all_pix2world(pixel_x, pixel_y, 1)
        depth_grid = np.full((len(y_edges) - 1, len(x_edges) - 1), np.nan, dtype=float)
        invalid_grid = np.zeros_like(depth_grid, dtype=bool)
        for row in band_rows:
            ix = int(np.searchsorted(x_edges, float(row['x_min']), side='left'))
            iy = int(np.searchsorted(y_edges, float(row['y_min']), side='left'))
            if ix >= depth_grid.shape[1] or iy >= depth_grid.shape[0]:
                continue
            status = str(row['status']).strip()
            depth = float(row['depth_5sigma_mag'])
            if status == 'ok' and np.isfinite(depth):
                depth_grid[iy, ix] = depth
            else:
                invalid_grid[iy, ix] = True
        combined_depth_maps[(cosmos, band)] = {
            'wcs': band_wcs,
            'exptime_s': exposure_seconds,
            'x_edges': x_edges + 0.5,
            'y_edges': y_edges + 0.5,
            'ra_edges': np.asarray(ra_edges, dtype=float),
            'dec_edges': np.asarray(dec_edges, dtype=float),
            'depth': depth_grid,
            'invalid': invalid_grid,
        }

COMBINED_DEPTH_RA_MIN = float(np.nanmin(all_ra))
COMBINED_DEPTH_RA_MAX = float(np.nanmax(all_ra))
COMBINED_DEPTH_DEC_MIN = float(np.nanmin(all_dec))
COMBINED_DEPTH_DEC_MAX = float(np.nanmax(all_dec))
if background_samples:
    background_values = np.concatenate(background_samples)
    background_vmin, background_vmax = np.nanpercentile(background_values, [1, 99])
else:
    background_vmin, background_vmax = 0.0, 1.0
COMBINED_DEPTH_BACKGROUND_NORM = Normalize(background_vmin, background_vmax, clip=True)

for band in BANDS_7DT:
    fields = sorted(combined_depth_available_by_band.get(band, set()))
    print(f'{band}: available COSMOS fields={fields or "none"}')
print(f'Combined RA range: [{COMBINED_DEPTH_RA_MIN:.5f}, {COMBINED_DEPTH_RA_MAX:.5f}] deg')
print(f'Combined Dec range: [{COMBINED_DEPTH_DEC_MIN:.5f}, {COMBINED_DEPTH_DEC_MAX:.5f}] deg')


# Use the same all-field/all-band magnitude limits as Proposal.ipynb.
# Round outward to half magnitudes; avoid per-band percentile clipping so
# the same AB limiting magnitude maps to the same color in every image.
COMBINED_DEPTH_ALL_VALUES = np.concatenate([
    entry['depth'][np.isfinite(entry['depth'])]
    for entry in combined_depth_maps.values()
])
if not COMBINED_DEPTH_ALL_VALUES.size:
    raise ValueError('No finite depth measurements for the shared color scale')
COMBINED_DEPTH_VMIN = float(np.floor(COMBINED_DEPTH_ALL_VALUES.min() * 2) / 2)
COMBINED_DEPTH_VMAX = float(np.ceil(COMBINED_DEPTH_ALL_VALUES.max() * 2) / 2)
if COMBINED_DEPTH_VMAX == COMBINED_DEPTH_VMIN:
    COMBINED_DEPTH_VMAX += 0.5
COMBINED_DEPTH_NORM = Normalize(COMBINED_DEPTH_VMIN, COMBINED_DEPTH_VMAX, clip=True)
COMBINED_DEPTH_CMAP = plt.get_cmap('viridis')
print(f'Shared depth color scale: {COMBINED_DEPTH_VMIN:.1f}–{COMBINED_DEPTH_VMAX:.1f} AB mag')


from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as patheffects

ra_edges = np.linspace(COMBINED_DEPTH_RA_MIN, COMBINED_DEPTH_RA_MAX, COMBINED_DEPTH_MAP_NX + 1)
dec_edges = np.linspace(COMBINED_DEPTH_DEC_MIN, COMBINED_DEPTH_DEC_MAX, COMBINED_DEPTH_MAP_NY + 1)
ra_centers = 0.5 * (ra_edges[:-1] + ra_edges[1:])
dec_centers = 0.5 * (dec_edges[:-1] + dec_edges[1:])
world_ra, world_dec = np.meshgrid(ra_centers, dec_centers)

detection_coverage = np.zeros(world_ra.shape, dtype=bool)
for cosmos in COMBINED_DEPTH_COSMOS:
    footprint = combined_depth_footprints[cosmos]
    pixel_x, pixel_y = footprint['wcs'].all_world2pix(world_ra, world_dec, 1)
    detection_coverage |= (pixel_x >= 0.5) & (pixel_x < footprint['nx_used'] + 0.5) & (pixel_y >= 0.5) & (pixel_y < footprint['ny_used'] + 0.5)


def proposal_draw_depth(fig, ax, band):
    """Draw the current depth-analysis recipe on one of the figure's axes.

    Use native axes with the shared all-band normalization, display-only Field
    numbering, and the original field/band coadd's FITS EXPTIME annotation.
    """
    available_fields = sorted(combined_depth_available_by_band.get(band, set()))
    selected_depth = np.full(world_ra.shape, np.nan, dtype=float)
    selected_field = np.zeros(world_ra.shape, dtype=np.int8)
    band_coverage = np.zeros(world_ra.shape, dtype=bool)
    for cosmos in available_fields:
        depth_map = combined_depth_maps[(cosmos, band)]
        pixel_x, pixel_y = depth_map['wcs'].all_world2pix(world_ra, world_dec, 1)
        ix = np.searchsorted(depth_map['x_edges'], pixel_x, side='right') - 1
        iy = np.searchsorted(depth_map['y_edges'], pixel_y, side='right') - 1
        inside = (ix >= 0) & (ix < depth_map['depth'].shape[1]) & (iy >= 0) & (iy < depth_map['depth'].shape[0])
        band_coverage |= inside
        candidate = np.full(world_ra.shape, np.nan, dtype=float)
        candidate[inside] = depth_map['depth'][iy[inside], ix[inside]]
        use_candidate = np.isfinite(candidate) & (~np.isfinite(selected_depth) | (candidate > selected_depth))
        selected_depth[use_candidate] = candidate[use_candidate]
        selected_field[use_candidate] = cosmos
    finite_depths = selected_depth[np.isfinite(selected_depth)]
    depth_vmin, depth_vmax = COMBINED_DEPTH_VMIN, COMBINED_DEPTH_VMAX
    depth_norm, depth_cmap = COMBINED_DEPTH_NORM, COMBINED_DEPTH_CMAP

    # Match both depth panels with a white background and black annotations.
    # Keep viridis and its AB-magnitude normalization unchanged scientifically.
    light_panel = band in ('m425', 'm650')
    text_color = 'black' if light_panel else 'white'
    box_color = 'white' if light_panel else 'black'
    ax.set_facecolor('white' if light_panel else '#252525')
    for cosmos in COMBINED_DEPTH_COSMOS:
        footprint = combined_depth_footprints[cosmos]
        ax.imshow(
            footprint['background'], origin='lower', cmap=('gray_r' if light_panel else 'gray'), norm=COMBINED_DEPTH_BACKGROUND_NORM,
            extent=[footprint['ra_max'], footprint['ra_min'], footprint['dec_min'], footprint['dec_max']],
            interpolation='nearest', alpha=COMBINED_DEPTH_BACKGROUND_ALPHA, aspect='auto', zorder=0,
        )

    ax.pcolormesh(
        ra_edges, dec_edges, np.ma.masked_invalid(selected_depth), cmap=depth_cmap,
        # Opaque colors match Proposal and do not mix with the grayscale background.
        norm=depth_norm, shading='flat', alpha=1.0, zorder=2, rasterized=True,
    )
    unavailable_only = detection_coverage & ~band_coverage
    insufficient_only = band_coverage & ~np.isfinite(selected_depth)
    if np.any(unavailable_only):
        unavailable_hatch = ax.contourf(
            world_ra, world_dec, unavailable_only.astype(float), levels=[0.5, 1.5],
            colors=('white' if light_panel else 'none'), hatches=['///'],
            alpha=(1.0 if light_panel else 0.55), zorder=3)
        if light_panel:
            # Explicit opaque white avoids alpha compositing the original
            # detection mosaic into a gray missing-band region in PDF viewers.
            unavailable_hatch.set_edgecolor('0.65')
            unavailable_hatch.set_linewidth(0)
    if np.any(insufficient_only):
        insufficient_hatch = ax.contourf(
            world_ra, world_dec, insufficient_only.astype(float), levels=[0.5, 1.5],
            colors=('white' if light_panel else 'none'), hatches=['..'],
            alpha=(1.0 if light_panel else 0.42), zorder=4)
        if light_panel:
            # Dots retain their original meaning: no valid depth in this cell.
            insufficient_hatch.set_edgecolor('0.65')
            insufficient_hatch.set_linewidth(0)

    for cosmos in COMBINED_DEPTH_COSMOS:
        footprint = combined_depth_footprints[cosmos]
        is_available = cosmos in available_fields
        field_number = COMBINED_DEPTH_DISPLAY_FIELD[cosmos]
        if is_available:
            # Exposure is specific to this original field AND this band.
            # No averaging across overlapping fields or filename-based values.
            exposure_seconds = combined_depth_maps[(cosmos, band)]['exptime_s']
            label = f'Field {field_number}\n({exposure_seconds:g}s)'
        else:
            # A missing band has no EXPTIME; do not imply a zero-second image.
            label = f'Field {field_number}\n(unavailable)'
        ax.text(
            np.mean(footprint['ra_corners']), np.mean(footprint['dec_corners']), label,
            ha='center', va='center', color=text_color, alpha=0.96, fontsize=22, fontweight='bold',
            bbox=dict(facecolor=box_color, edgecolor='none', alpha=COMBINED_DEPTH_LABEL_BBOX_ALPHA, pad=2.5),
            path_effects=[patheffects.withStroke(linewidth=1.5, foreground=box_color, alpha=0.55)],
            zorder=6,
        )

    # Place the band name at the upper-left corner of the panel, using the
    # same white type and translucent black box as the field annotations.
    ax.text(0.04, 0.96, band, transform=ax.transAxes,
            ha='left', va='top', color=text_color, alpha=0.96,
            fontsize=22, fontweight='bold',
            bbox=dict(facecolor=box_color, edgecolor='none',
                      alpha=COMBINED_DEPTH_LABEL_BBOX_ALPHA, pad=2.5),
            path_effects=[patheffects.withStroke(linewidth=1.5, foreground=box_color, alpha=0.55)],
            zorder=7)

    depth_mappable = ScalarMappable(norm=depth_norm, cmap=depth_cmap)
    depth_mappable.set_array(finite_depths if finite_depths.size else np.array([depth_vmin, depth_vmax]))
    colorbar = fig.colorbar(depth_mappable, ax=ax, pad=0.02)
    colorbar.set_label(r'$5\sigma$ limiting magnitude [AB mag]', fontsize=22.4)
    colorbar.ax.tick_params(labelsize=18)
    ax.set_xlim(COMBINED_DEPTH_RA_MAX + COMBINED_DEPTH_PAD_DEG, COMBINED_DEPTH_RA_MIN - COMBINED_DEPTH_PAD_DEG)
    ax.set_ylim(COMBINED_DEPTH_DEC_MIN - COMBINED_DEPTH_PAD_DEG, COMBINED_DEPTH_DEC_MAX + COMBINED_DEPTH_PAD_DEG)
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('Dec [deg]')
    ax.grid(alpha=0.16, linewidth=0.6)

# Assemble the new single figure; its main axes are exactly one row by three.
PROPOSAL_THREE_OUTPUT = FINAL_OUTPUT
PROPOSAL_THREE_OUTPUT.mkdir(parents=True, exist_ok=True)
with plt.rc_context():
    scl_plt.default()
    # Embed TrueType fonts; axes, curves, hatches, and text remain vector.
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(1, 3, figsize=(25, 7), layout='constrained')
    hist_ax = axes[0]
    # Match the latest notebook 2.3 cumulative figure: bin-center connected
    # lines and circular band markers, rather than step/stair rendering.
    hist_ax.plot(CUM_BIN_CENTERS, CUM_REFERENCE_COUNTS[1:],
                 color='black', linestyle='--', linewidth=1.8,
                 label='Subaru/HSC', zorder=4)
    for band, wavelength in zip(CUM_BANDS, CUM_WAVELENGTH):
        hist_ax.plot(CUM_BIN_CENTERS, CUM_BAND_COUNTS[band][1:],
                     linewidth=1.25,
                     color=CUM_CMAP(CUM_NORM(wavelength)))
    hist_ax.set_xlabel(r'$r$', fontsize=22)
    hist_ax.set_ylabel(r'$N_{\mathrm{source}}$', fontsize=22)
    hist_ax.tick_params(axis='both', labelsize=18)
    hist_ax.set_xlim(17, 21)
    visible = (CUM_BIN_CENTERS >= 17) & (CUM_BIN_CENTERS <= 21)
    boundary = np.interp([17, 21], CUM_BIN_CENTERS, CUM_REFERENCE_COUNTS[1:])
    visible_max = max(1, np.max(CUM_REFERENCE_COUNTS[1:][visible]), np.max(boundary))
    hist_ax.set_ylim(0, visible_max * 1.08)
    hist_ax.grid(alpha=0.25)
    hist_ax.legend(loc='upper left', fontsize=30)
    mappable = plt.cm.ScalarMappable(norm=CUM_NORM, cmap=CUM_CMAP)
    colorbar = fig.colorbar(mappable, ax=hist_ax, orientation='vertical', pad=0.02, fraction=0.05)
    colorbar.set_label('7DT central wavelength [nm]', fontsize=22.4)
    colorbar.ax.tick_params(labelsize=18)
    for ax, band in zip(axes[1:], ('m425', 'm650')):
        proposal_draw_depth(fig, ax, band)
    # Enlarge only the coordinate/magnitude axis labels; field and band
    # annotations are independent text artists and retain their current size.
    # Match each panel's coordinate tick numbers to its adjacent colorbar.
    hist_ax.tick_params(axis='both', labelsize=18)
    for ax in axes[1:]:
        ax.tick_params(axis='both', labelsize=18)
    for ax in axes:
        ax.set_xlabel(ax.get_xlabel(), fontsize=ax.xaxis.label.get_fontsize() * 1.7)
        ax.set_ylabel(ax.get_ylabel(), fontsize=ax.yaxis.label.get_fontsize() * 1.7)
    stem = PROPOSAL_THREE_OUTPUT / 'figure1_baseline_flagcut_depth_white_600dpi'
    fig.savefig(stem.with_suffix('.png'), dpi=600, bbox_inches='tight')
    fig.savefig(stem.with_suffix('.pdf'), dpi=600, bbox_inches='tight')
    # The preview is separately rendered from the PDF during visual QA.
    plt.close(fig)
print(f'Saved {stem.with_suffix(".png")}')
print(f'Saved {stem.with_suffix(".pdf")}')


# Preserve the exact scientific input identities alongside the deliverables.
def sha256(path):
    """Hash the small catalog/table inputs in chunks without altering them."""
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

inputs = [CUM_PATH, CUM_ROOT / 'COSMOS_data/COSMOS_finalcat.fits']
inputs += [_combined_depth_csv_path(field) for field in COMBINED_DEPTH_COSMOS]
inputs += [CUM_ROOT / 'output/experiments/aperture_flag_ablation/baseline/Notebooks/Proposal.ipynb']
manifest = {
    'selection': 'Archived pre-ablation results, retaining aper3_flags == 0',
    'recipe': 'Archived Proposal.ipynb cells 37-39, with display/output changes only',
    'png_dpi': 600,
    'pdf': 'Vector curves, axes, hatches and annotations; 600-dpi raster depth layer; embedded TrueType fonts',
    'display_change': 'Both depth panels use light backgrounds and black annotations; viridis unchanged',
    'counts_verified': 'Exact equality with all columns of archived cumulative-count CSV',
    'input_sha256': {str(path.relative_to(CUM_ROOT)): sha256(path) for path in inputs},
    'outputs': [str(stem.with_suffix(suffix).relative_to(CUM_ROOT)) for suffix in ('.png', '.pdf')],
}
(FINAL_OUTPUT / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
(FINAL_OUTPUT / 'caption_note.txt').write_text(
    'Caption revision for the right panel: Replace "The black area indicates there is no image in there." '
    'with "The hatched area indicates that no image is available in this band."\n'
    'Hatched areas indicate an unavailable band; dotted areas indicate cells without a valid depth estimate.\n')
