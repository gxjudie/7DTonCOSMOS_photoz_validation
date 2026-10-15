"""Shared two-page A4 source inspection for notebooks 2.2 and 2.3.

The plotting calculations originate in notebook 2.2. Catalog adaptation is limited
by field: global IDs remain global, while positions, photometry, and PSFs must all
refer to the same COSMOS image. Display and PDF export use identical figures.
"""
from functools import lru_cache
from pathlib import Path
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.table import Table as FigureTable
from matplotlib.textpath import TextPath
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

BANDS_7DT = tuple(f'm{i}' for i in range(400, 900, 25))
A4_WIDTH, A4_HEIGHT = 595.275590551, 841.88976378
NAVY, MUTED, PALE = '#172b4d', '#546274', '#f1f5f9'


class SourceInfoResult(dict):
    """Return tables and report figures without dumping arrays in notebook output."""
    def __repr__(self):
        return f"SourceInfoResult(master_id={self['master_id']}, pdf_path={self.get('pdf_path')})"

    def _repr_pretty_(self, printer, cycle):
        printer.text(repr(self))


def get_image_path(project_root, cosmos, band):
    """Require a unique science image rather than silently choosing a file."""
    paths = sorted((project_root / f'COSMOS_data/COSMOS{cosmos}').glob(f'*{band}*.com.fits'))
    if len(paths) != 1:
        raise FileNotFoundError(f'Expected one COSMOS{cosmos} {band} image; found {len(paths)}')
    return paths[0]


@lru_cache(maxsize=1)
def _field_images(project_root, cosmos):
    """Cache only the current field to avoid retaining four sets of large images."""
    root = Path(project_root)
    data_dir = root / f'COSMOS_data/COSMOS{cosmos}'
    # Availability follows the existing notebook requirement: science plus weight.
    active = tuple(b for b in BANDS_7DT
                   if len(list(data_dir.glob(f'*{b}*.com.fits'))) == 1
                   and len(list(data_dir.glob(f'*{b}*.weight.fits'))) == 1)
    if not active:
        raise FileNotFoundError(f'No science/weight band pairs for COSMOS{cosmos}')
    reference = min(active, key=lambda b: abs(int(b[1:]) - 600))
    candidates = [data_dir / f'COSMOS{cosmos}_detectionimage.fits',
                  root / f'output/COSMOS{cosmos}/COSMOS{cosmos}_detectionimage.fits']
    detection_path = next((p for p in candidates if p.exists()), None)
    if detection_path is None:
        detection_path = get_image_path(root, cosmos, reference)
    images = {b: fits.getdata(get_image_path(root, cosmos, b), memmap=True) for b in active}
    return active, images, fits.getdata(detection_path, memmap=True)


def catalog_for_field(catalog, field):
    """Expose field-local measurement names while retaining global master IDs.

    This is a column projection of the supplied persisted catalog, not a new join
    or a reload of the original field master. Reconciled SDSS stays common; field
    FINALCAT match distances/status take precedence over shared external columns.
    """
    prefix = f'7dt_cosmos{field}_'
    if f'{prefix}matched' not in catalog.colnames:
        return catalog
    subset = catalog[np.asarray(np.ma.filled(catalog[f'{prefix}matched'], False), dtype=bool)]
    view = Table(masked=True)
    for name in subset.colnames:
        if name.startswith(prefix) and name != f'{prefix}matched':
            view[name[len(prefix):]] = subset[name].copy()
    view['field_master_id'] = view['master_id'].copy()
    view['master_id'] = subset['master_id'].copy()
    for name in subset.colnames:
        if (name.startswith('finalcat_') or name.startswith('sdss_') or name == 'has_sdss') and name not in view.colnames:
            view[name] = subset[name].copy()
    return view


def _page_header(fig, title, source_id, cosmos, page):
    """Use physical A4 dimensions and consistent, compact page furniture."""
    fig.text(.034, .967, title, fontsize=11, color=NAVY)
    fig.text(.966, .97, f'master_id {source_id} | COSMOS{cosmos}',
             fontsize=6.5, ha='right', color=MUTED)
    fig.text(.034, .015, 'Source inspection', fontsize=5.5, color=MUTED)
    fig.text(.966, .015, f'{page} / 2', fontsize=5.5, ha='right', color=MUTED)


def _image_page(source_id, cosmos):
    """Create one cutout grid and six axes below it; SED spans the right column."""
    fig = plt.figure(figsize=(A4_WIDTH / 72, A4_HEIGHT / 72))
    _page_header(fig, 'Images and diagnostics', source_id, cosmos, 1)
    fig.text(.034, .938, '20-band cutouts', fontsize=7, color=NAVY)
    bands = fig.add_gridspec(4, 5, left=.065, right=.978, bottom=.448, top=.915,
                            wspace=.27, hspace=.28)
    band_axes = np.array([[fig.add_subplot(bands[r, c]) for c in range(5)] for r in range(4)])
    fig.text(.034, .424, 'Detection image', fontsize=7, color=NAVY)
    fig.text(.426, .424, 'Population diagnostics and SED', fontsize=7, color=NAVY)
    image_ax = fig.add_axes([.075, .165, .29, .24])
    diagnostics = fig.add_gridspec(2, 2, left=.455, right=.975, bottom=.185, top=.398,
                                  wspace=.43, hspace=.60)
    axes = np.empty((2, 3), dtype=object)
    axes[0, 0] = fig.add_subplot(diagnostics[0, 0])
    axes[0, 1] = fig.add_subplot(diagnostics[0, 1])
    axes[0, 2] = fig.add_subplot(diagnostics[1, 0])
    axes[1, 0] = fig.add_subplot(diagnostics[1, 1])
    # The wider SED is an intentional change from the approved compact prototype.
    axes[1, 1] = fig.add_axes([.455, .055, .52, .090])
    axes[1, 2] = None
    return fig, image_ax, band_axes, axes


def _style_image_page(fig, image_ax, band_axes, axes):
    """Scale labels and symbols to their A4 panel sizes, keeping plotted values."""
    for ax in band_axes.flat:
        ax.title.set_fontsize(4.7)
        ax.tick_params(labelsize=3.3, length=1.5, pad=1)
        for collection in ax.collections:
            if hasattr(collection, 'get_sizes'):
                collection.set_sizes(np.asarray(collection.get_sizes()) * .20)
    image_ax.title.set_fontsize(5.4)
    image_ax.xaxis.label.set_size(5)
    image_ax.yaxis.label.set_size(5)
    image_ax.tick_params(labelsize=4.2, length=2, pad=1)
    image_ax.legend(loc='upper right', fontsize=3.8, markerscale=.5, borderpad=.3)
    for ax in [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]]:
        ax.title.set_fontsize(5.1)
        ax.xaxis.label.set_size(4.8)
        ax.yaxis.label.set_size(4.8)
        ax.tick_params(which='both', labelsize=4, length=2, pad=1)
        ax.yaxis.get_offset_text().set_fontsize(4)
        for line in ax.lines:
            line.set_markersize(min(line.get_markersize(), 2.5))
            line.set_linewidth(min(line.get_linewidth(), .65))
        for collection in ax.collections:
            if hasattr(collection, 'get_sizes'):
                collection.set_sizes(np.asarray(collection.get_sizes()) * .30)
        if ax.get_legend() is not None:
            ax.legend(fontsize=3.2, markerscale=.5, borderpad=.25, handlelength=1.1,
                      labelspacing=.25)


def _format_value(raw):
    """Never round integer identifiers; match the notebook's six-digit display."""
    if np.ma.is_masked(raw):
        return '--'
    if isinstance(raw, (bytes, np.bytes_)):
        return raw.decode(errors='replace').strip()
    if isinstance(raw, (float, np.floating)):
        return f'{raw:.6g}' if np.isfinite(raw) else '--'
    return str(raw)


@lru_cache(maxsize=4096)
def _wrap_cell(text, width, fontsize):
    """Measure actual glyph widths and preserve deliberate multiline headers."""
    count = max(3, int((width-4) / (fontsize * .45)))
    while count >= 2:
        lines = []
        for paragraph in text.split('\n'):
            lines.extend(textwrap.wrap(paragraph, width=count, break_long_words=True,
                                       break_on_hyphens=False) or [''])
        if all(not line or TextPath((0, 0), line, size=fontsize).get_extents().width <= width-4
               for line in lines):
            return '\n'.join(lines)
        count -= 1
    return '\n'.join(text)


def _table_layout(rows, widths, fontsize, headers=None):
    """Wrap cells before computing heights so long flags cannot overlap rows."""
    data = ([] if headers is None else [headers]) + rows
    wrapped, heights = [], []
    for row in data:
        cells = [_wrap_cell(_format_value(raw), width, fontsize) for raw, width in zip(row, widths)]
        wrapped.append(cells)
        heights.append(max(cell.count('\n') + 1 for cell in cells) * fontsize * 1.2 + 4)
    return wrapped, heights


def _draw_table(fig, rows, widths, x, top, fontsize, headers=None):
    """Draw compact vector text tables without depending on a document library."""
    wrapped, heights = _table_layout(rows, widths, fontsize, headers)
    height, width = sum(heights), sum(widths)
    ax = fig.add_axes([x / A4_WIDTH, (top-height) / A4_HEIGHT,
                       width / A4_WIDTH, height / A4_HEIGHT])
    ax.axis('off')
    table = FigureTable(ax, bbox=[0, 0, 1, 1])
    for r, (row, rh) in enumerate(zip(wrapped, heights)):
        for col, (text, cw) in enumerate(zip(row, widths)):
            header = headers is not None and r == 0
            cell = table.add_cell(r, col, cw/width, rh/height, text=text, loc='left',
                                  facecolor=NAVY if header else (PALE if r % 2 == 0 else 'white'),
                                  edgecolor='#dce3ec')
            cell.set_linewidth(.2)
            cell.PAD = 2 / cw
            cell.get_text().set_color('white' if header else NAVY)
            cell.get_text().set_fontsize(fontsize)
    table.auto_set_font_size(False)
    ax.add_table(table)
    return top-height


def _table_page(summary, band_table, finalcat_table, sdss_table, source_id, cosmos):
    """Fit all scientific tables on the second portrait A4 page without omission."""
    fig = plt.figure(figsize=(A4_WIDTH / 72, A4_HEIGHT / 72))
    _page_header(fig, 'Source data', source_id, cosmos, 2)
    source_rows = [(name, _format_value(summary[name][0])) for name in summary.colnames]
    half = (len(source_rows)+1)//2
    summary_rows = [list(source_rows[i]) + (list(source_rows[i+half]) if i+half<len(source_rows) else ['', ''])
                    for i in range(half)]
    columns = ['band', 'wavelength_nm', 'detected', 'source_id', 'mag', 'mag_err', 'snr',
               'ellipticity', 'flux_radius_pix', 'derived_fwhm_pix', 'derived_fwhm_arcsec',
               'band_fwhm_pix', 'band_fwhm_arcsec', 'extendicity']
    headers = ['Band', 'Wave.\nnm', 'Detected', 'Source\nID', 'Mag', 'Mag\nerror', 'S/N',
               'Ellipticity', 'Flux radius\npix', 'Derived\nFWHM\npix', 'Derived\nFWHM\narcsec',
               'Band\nFWHM\npix', 'Band\nFWHM\narcsec', 'Extendicity']
    band_rows = [[_format_value(row[name]) for name in columns] for row in band_table]
    final_rows = [[_format_value(v) for v in row] for row in finalcat_table]
    sdss_rows = [[_format_value(v) for v in row] for row in sdss_table]
    width = A4_WIDTH-40
    source_widths = [150, 65, 235, width-450]
    weights = np.array([26,30,30,34,36,42,36,39,40,41,41,41,41,42], dtype=float)
    band_widths = list(weights / weights.sum() * width)
    left_width = 251
    right_x, right_width = 20+left_width+12, width-left_width-12
    # Rare long metadata/flags reduce all table fonts together, never drop rows.
    for scale in np.arange(1., .69, -.025):
        summary_height = sum(_table_layout(summary_rows, source_widths, 6.8*scale)[1])
        band_height = sum(_table_layout(band_rows, band_widths, 6.2*scale, headers)[1])
        external_height = max(sum(_table_layout(final_rows, [135,left_width-135], 6.8*scale)[1]),
                              sum(_table_layout(sdss_rows, [100,right_width-100], 6.8*scale)[1]))
        if A4_HEIGHT-58-summary_height-24-band_height-24-external_height >= 42:
            break
    else:
        plt.close(fig)
        raise ValueError('Source metadata is too long for a legible two-page A4 report.')
    def section(title, top, x=20):
        fig.text(x/A4_WIDTH, top/A4_HEIGHT, title, fontsize=7.5, color=NAVY)
        return top-7
    top = section('Source information and morphology summary', A4_HEIGHT-51)
    bottom = _draw_table(fig, summary_rows, source_widths, 20, top, 6.8*scale)
    top = section('7DT photometry and morphology', bottom-17)
    bottom = _draw_table(fig, band_rows, band_widths, 20, top, 6.2*scale, headers)
    top = bottom-17
    _draw_table(fig, final_rows, [135,left_width-135], 20, section('FINALCAT information', top), 6.8*scale)
    _draw_table(fig, sdss_rows, [100,right_width-100], right_x, section('SDSS information', top, right_x), 6.8*scale)
    fig.text(.034, .034, '--: missing / invalid. Derived FWHM = 2 x flux radius; extendicity = derived FWHM / band FWHM. Pixels: FITS 1-based.',
             fontsize=5, color=MUTED)
    return fig


def plot_source_info(
    cat, image=None, pixel_range=40, source_id=None, *,
    info_catalog=None, population_catalog=None, diagnostic_band='m600', show=True,
    save_pdf=False, output_dir=None, project_root=None, field=None,
    active_bands=None, band_images=None, detection_image=None,
):
    """Display a two-page A4 source report and optionally save the same pages.

    Existing catalog selection and scientific calculations follow notebook 2.2.
    ``save_pdf=True`` chooses an ID-based filename; a path chooses an explicit PDF.
    ``show=False`` suppresses notebook display and closes report figures after save.
    Legacy figure keys alias page one; use ``report_figures`` for both report pages.
    """
    project_root = Path(project_root or Path.cwd()).resolve()
    info_catalog = cat if info_catalog is None else info_catalog
    population_catalog = cat if population_catalog is None else population_catalog
    if 'master_id' not in cat.colnames or len(cat) == 0:
        raise ValueError('The input catalog must contain sources and master_id.')
    ids = np.asarray(cat['master_id'], dtype=np.int64)
    source_id = int(ids[np.random.default_rng().integers(len(ids))]) if source_id is None else int(source_id)
    selected = np.flatnonzero(ids == source_id)
    if len(selected) != 1:
        raise KeyError(f'master_id={source_id} must occur exactly once in cat; found {len(selected)}.')
    original_row = cat[int(selected[0])]
    unified = 'selected_COSMOS' in cat.colnames and '7dt_cosmos1_matched' in cat.colnames
    if unified:
        # Keep the global identity while using only this field's pixels and photometry.
        field = int(original_row['selected_COSMOS'])
        active_bands = band_images = detection_image = None
        cat = catalog_for_field(cat, field)
        info_catalog = catalog_for_field(info_catalog, field)
        population_catalog = catalog_for_field(population_catalog, field)
    COSMOS = int(1 if field is None else field)
    needs_images = band_images is None or (image is None and detection_image is None)
    loaded_bands, loaded_images, loaded_detection = (
        _field_images(str(project_root), COSMOS) if needs_images else (None, None, None)
    )
    ACTIVE_BANDS = tuple(active_bands or loaded_bands or band_images.keys())
    band_images = loaded_images if band_images is None else band_images
    detection_image = loaded_detection if detection_image is None else detection_image
    image = detection_image if image is None else image
    N_ACTIVE_BANDS = len(ACTIVE_BANDS)
    REFERENCE_BAND = min(ACTIVE_BANDS, key=lambda b: abs(int(b[1:]) - 600))
    MATCH_DIR = project_root / f'output/COSMOS{COSMOS}/singlephot/diagnostics'
    if diagnostic_band not in ACTIVE_BANDS:
        # Automatic fallback applies only to the usual m600 default in fields without m600.
        if diagnostic_band == 'm600':
            diagnostic_band = REFERENCE_BAND
        else:
            raise ValueError(f'diagnostic_band must be one of {ACTIVE_BANDS}.')
    if len(cat) == 0:
        raise ValueError('The input catalog is empty.')
    required = {'master_id', 'master_x', 'master_y'}
    missing = required.difference(cat.colnames)
    if missing:
        raise KeyError(f'Catalog is missing required columns: {sorted(missing)}')
    if pixel_range is None or not np.isfinite(pixel_range) or pixel_range <= 0:
        raise ValueError('pixel_range must be a positive number.')
    info_catalog = cat if info_catalog is None else info_catalog
    population_catalog = cat if population_catalog is None else population_catalog
    info_ids = np.asarray(info_catalog['master_id'], dtype=np.int64)
    cat_ids = np.asarray(cat['master_id'], dtype=np.int64)
    if source_id is None:
        selected_id = int(cat_ids[np.random.default_rng().integers(len(cat))])
    else:
        selected_id = int(source_id)
    matches = np.flatnonzero(info_ids == selected_id)
    cat_matches = np.flatnonzero(cat_ids == selected_id)
    if len(matches) != 1:
        raise KeyError(f'master_id={selected_id} must occur exactly once in info_catalog; found {len(matches)}.')
    if len(cat_matches) != 1:
        raise KeyError(f'master_id={selected_id} must occur exactly once in cat; found {len(cat_matches)}.')
    row = info_catalog[int(matches[0])]
    image_row = cat[int(cat_matches[0])]

    def value(name, default=np.nan):
        if name not in row.colnames:
            return default
        raw = row[name]
        if np.ma.is_masked(raw):
            return default
        return raw.item() if hasattr(raw, 'item') else raw

    def number(name):
        try:
            return float(value(name))
        except (TypeError, ValueError):
            return np.nan

    def display_value(name):
        raw = value(name, np.ma.masked)
        if np.ma.is_masked(raw):
            return '--'
        if isinstance(raw, (bytes, np.bytes_, str)):
            text = raw.decode(errors='replace').strip() if isinstance(raw, (bytes, np.bytes_)) else str(raw).strip()
            return '--' if text in {'', 'NN', 'nan'} else text
        try:
            numeric = float(raw)
        except (TypeError, ValueError):
            return str(raw)
        if not np.isfinite(numeric) or numeric in {-9.0, -99.0, -999.0}:
            return '--'
        return str(int(raw)) if isinstance(raw, (int, np.integer)) else f'{numeric:.6g}'

    selected_x = float(image_row['master_x'])
    selected_y = float(image_row['master_y'])
    if not np.isfinite(selected_x) or not np.isfinite(selected_y):
        raise ValueError('The selected source has a non-finite pixel position in cat.')
    summary = Table()
    summary['master_id'] = [selected_id]
    summary['master_class'] = [str(value('master_class', ''))]
    summary['master_ra'] = [number('master_ra')]
    summary['master_dec'] = [number('master_dec')]
    summary['master_x'] = [selected_x]
    summary['master_y'] = [selected_y]
    summary['n_7dt_bands'] = [int(value('n_7dt_bands', 0))]
    summary['has_finalcat'] = [bool(value('has_finalcat', False))]
    summary['has_sdss'] = [bool(value('has_sdss', False))]
    report_fig, image_ax, band_axes, axes = _image_page(source_id, COSMOS)
    image_fig, image_ax, selected_id = plot_source_image(
        cat, image=image, pixel_range=pixel_range, source_id=selected_id, ax=image_ax
    )

    x = int(round(selected_x - 1))
    y = int(round(selected_y - 1))
    samples = []
    for band in ACTIVE_BANDS:
        flat = band_images[band].reshape(-1)
        step = max(1, len(flat) // 5000)
        sample = np.asarray(flat[::step][:5000], dtype=float)
        samples.append(sample[np.isfinite(sample)])
    pooled = np.concatenate(samples)
    if pooled.size == 0:
        raise ValueError('The 7DT images contain no finite pixels for ZScale.')
    vmin, vmax = ZScaleInterval(n_samples=len(pooled)).get_limits(pooled)
    band_norm = ImageNormalize(vmin=vmin, vmax=vmax)
    xs = np.ma.filled(cat['master_x'], np.nan).astype(float)
    ys = np.ma.filled(cat['master_y'], np.nan).astype(float)
    classes = np.asarray(cat['master_class']).astype(str) if 'master_class' in cat.colnames else np.full(len(cat), 'source')
    class_colors = {'matched': 'purple', '7dt_only': 'tab:red', 'finalcat_only': 'tab:blue'}
    band_fig = report_fig
    for band, ax in zip(BANDS_7DT, band_axes.flat):
        if band not in ACTIVE_BANDS:
            ax.set_title(f'{band} (not observed)', fontsize=5)
            ax.text(.5, .5, 'Not observed', transform=ax.transAxes, ha='center', va='center', fontsize=4)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        band_image = band_images[band]
        half = int(round(pixel_range))
        xmin, xmax = max(0, x - half), min(band_image.shape[1], x + half + 1)
        ymin, ymax = max(0, y - half), min(band_image.shape[0], y + half + 1)
        crop = band_image[ymin:ymax, xmin:xmax]
        ax.imshow(crop, origin='lower', cmap='gray', norm=band_norm,
                  extent=[xmin + 0.5, xmax + 0.5, ymin + 0.5, ymax + 0.5])
        visible = np.isfinite(xs) & np.isfinite(ys) & (xs - 1 >= xmin) & (xs - 1 < xmax) & (ys - 1 >= ymin) & (ys - 1 < ymax)
        for source_class in np.unique(classes[visible]):
            selected = visible & (classes == source_class)
            ax.scatter(xs[selected], ys[selected], s=10, facecolors='none', edgecolors=class_colors.get(source_class, 'tab:red'), linewidths=0.4)
        detected = f'{band}_source_id' in row.colnames and not np.ma.is_masked(row[f'{band}_source_id']) and int(row[f'{band}_source_id']) >= 0
        ax.scatter(selected_x, selected_y, marker='+', s=80, color='cyan' if detected else 'gray')
        ax.set_title(f'{band} ({"detected" if detected else "not detected"})', fontsize=9)
        ax.tick_params(labelsize=7)
    uncertainty_path = MATCH_DIR / f'cosmos{COSMOS}_7dt_singleband_band_uncertainty.fits'
    band_uncertainty = Table.read(uncertainty_path) if uncertainty_path.exists() else None
    reference_fwhm = {} if band_uncertainty is None else {
        (band.decode().strip() if isinstance(band, (bytes, np.bytes_)) else str(band).strip()): float(fwhm)
        for band, fwhm in zip(band_uncertainty['band'], band_uncertainty['fwhm_pix'])
    }
    pixel_scale = float(proj_plane_pixel_scales(WCS(fits.getheader(get_image_path(project_root, COSMOS, REFERENCE_BAND))))[0] * 3600.0)
    band_rows = []
    wavelengths, mags, magerrs, snrs, extendicities = [], [], [], [], []
    for band in ACTIVE_BANDS:
        source_raw = value(f'{band}_source_id', np.ma.masked)
        detected = not np.ma.is_masked(source_raw) and int(source_raw) >= 0
        source_band_id = int(source_raw) if detected else -1
        mag, err = number(f'{band}_mag'), number(f'{band}_mag_err')
        mag = mag if detected and np.isfinite(mag) and -50 < mag < 50 else np.nan
        err = err if detected and np.isfinite(err) and err > 0 else np.nan
        # Small-error magnitude propagation: sigma_mag = (2.5 / ln(10)) / SNR.
        snr = 1.0857 / err if np.isfinite(err) else np.nan
        radius = number(f'{band}_flux_radius') if detected else np.nan
        radius = radius if np.isfinite(radius) and radius > 0 else np.nan
        ellipticity = number(f'{band}_ellipticity') if detected else np.nan
        ellipticity = ellipticity if np.isfinite(ellipticity) and 0 <= ellipticity <= 1 else np.nan
        # Preserve the original size proxy; this is not a fitted source FWHM.
        source_fwhm = 2.0 * radius if np.isfinite(radius) and radius > 0 else np.nan
        psf_fwhm = reference_fwhm.get(band, np.nan)
        extendicity = source_fwhm / psf_fwhm if np.isfinite(source_fwhm) and np.isfinite(psf_fwhm) and psf_fwhm > 0 else np.nan
        if np.isfinite(extendicity): extendicities.append(extendicity)
        band_rows.append((band, int(band[1:]), detected, source_band_id, mag, err, snr, ellipticity, radius, source_fwhm, source_fwhm * pixel_scale, psf_fwhm, psf_fwhm * pixel_scale, extendicity))
        wavelengths.append(int(band[1:])); mags.append(mag); magerrs.append(err); snrs.append(snr)
    band_table = Table(rows=band_rows, names=['band', 'wavelength_nm', 'detected', 'source_id', 'mag', 'mag_err', 'snr', 'ellipticity', 'flux_radius_pix', 'derived_fwhm_pix', 'derived_fwhm_arcsec', 'band_fwhm_pix', 'band_fwhm_arcsec', 'extendicity'])
    summary['extendicity_n_bands'] = [len(extendicities)]
    summary['extendicity_median'] = [float(np.nanmedian(extendicities)) if extendicities else np.nan]
    summary['extendicity_mad'] = [float(np.nanmedian(np.abs(np.asarray(extendicities) - np.nanmedian(extendicities)))) if extendicities else np.nan]

    def external_table(flag_name, fields, status_value):
        if not bool(value(flag_name, False)):
            diagnostic_names = {'has_sdss', 'sdss_match_status', 'sdss_match_method', 'sdss_n_candidates', 'sdss_n_master_candidates', 'finalcat_match_status'}
            rows = [(label, display_value(name)) for name, label in fields if name in diagnostic_names and name in row.colnames]
            rows.append(('status', status_value))
            return Table(rows=rows, names=['field', 'value'])
        rows = [(label, display_value(name)) for name, label in fields if name in row.colnames]
        return Table(rows=rows, names=['field', 'value'])
    finalcat_table = external_table('has_finalcat', [
        ('finalcat_match_status', 'match status'), ('finalcat_ID', 'ID'), ('finalcat_sep_arcsec', 'separation arcsec'),
        ('finalcat_STELLARITY', 'STELLARITY'), ('finalcat_CREFF', 'circular effective radius'),
        ('finalcat_REFF', 'effective radius'), ('finalcat_M13_STAR', 'M13 star flag'), ('finalcat_M13_SG', 'M13 star/galaxy flag'),
        ('finalcat_Z_TOT_Z', 'total redshift'), ('finalcat_Z_TOT_ZERR', 'redshift error'), ('finalcat_Z_TOT_ZSOURCE', 'redshift source'),
    ], 'No FINALCAT counterpart')
    sdss_table = external_table('has_sdss', [
        ('has_sdss', 'has SDSS'), ('sdss_match_status', 'match status'), ('sdss_match_method', 'match method'), ('sdss_objid', 'objID'),
        ('sdss_ra', 'RA deg'), ('sdss_dec', 'DEC deg'), ('sdss_sep_arcsec', 'separation arcsec'),
        ('sdss_n_candidates', 'SDSS candidates'), ('sdss_n_master_candidates', 'master candidates'),
        ('sdss_type', 'type'), ('sdss_clean', 'clean'), ('sdss_mode', 'mode'),
        ('sdss_flags', 'photo flags'), ('sdss_flags_text', 'photo flags text'),
        ('sdss_petroMag_u', 'Petrosian u'), ('sdss_petroMag_g', 'Petrosian g'), ('sdss_petroMag_r', 'Petrosian r'), ('sdss_petroMag_i', 'Petrosian i'), ('sdss_petroMag_z', 'Petrosian z'),
        ('sdss_petroMagErr_u', 'Petrosian u error'), ('sdss_petroMagErr_g', 'Petrosian g error'), ('sdss_petroMagErr_r', 'Petrosian r error'),
        ('sdss_petroMagErr_i', 'Petrosian i error'), ('sdss_petroMagErr_z', 'Petrosian z error'),
    ], 'No accepted SDSS counterpart (unmatched or ambiguous)')
    def phot(cat_, band, suffix='mag'):
        name = f'{band}_{suffix}'
        values = np.ma.filled(cat_[name], np.nan).astype(float) if name in cat_.colnames else np.full(len(cat_), np.nan)
        if suffix == 'mag': values[~np.isfinite(values) | (values <= -50) | (values >= 50)] = np.nan
        if suffix == 'mag_err': values[~np.isfinite(values) | (values <= 0)] = np.nan
        return values
    def source_phot(band, suffix='mag'):
        result = number(f'{band}_{suffix}')
        if suffix == 'mag' and (not np.isfinite(result) or result <= -50 or result >= 50): return np.nan
        if suffix == 'mag_err' and (not np.isfinite(result) or result <= 0): return np.nan
        return result
    if diagnostic_band not in ACTIVE_BANDS:
        raise ValueError(f'diagnostic_band must be one of {ACTIVE_BANDS}.')
    def selected_mark(ax, xx, yy, xerr=None, yerr=None):
        if np.isfinite(xx) and np.isfinite(yy):
            ax.errorbar(xx, yy, xerr=xerr, yerr=yerr, fmt='o', ms=8, color='red', markeredgecolor='black', label='selected source', zorder=5)

    m475, m600, m750 = [phot(population_catalog, b) for b in ('m475', 'm600', 'm750')]
    pop_fig = report_fig
    valid = np.isfinite(m475) & np.isfinite(m600) & np.isfinite(m750)
    axes[0, 0].scatter(m475[valid] - m750[valid], m600[valid], s=6, alpha=0.18, color='black', rasterized=True)
    e475, e600, e750 = [source_phot(band, 'mag_err') for band in ('m475', 'm600', 'm750')]
    cmd_xerr = np.hypot(e475, e750) if np.isfinite(e475) and np.isfinite(e750) else None
    selected_mark(axes[0, 0], source_phot('m475') - source_phot('m750'), source_phot('m600'), cmd_xerr, e600 if np.isfinite(e600) else None)
    axes[0, 0].set(xlabel='m475 - m750', ylabel='m600', title='Color-magnitude diagram'); axes[0, 0].invert_yaxis()
    if axes[0, 0].get_legend_handles_labels()[0]: axes[0, 0].legend(fontsize=8)
    axes[0, 1].scatter(m475[valid] - m600[valid], m600[valid] - m750[valid], s=6, alpha=0.18, color='black', rasterized=True)
    ccd_xerr = np.hypot(e475, e600) if np.isfinite(e475) and np.isfinite(e600) else None
    ccd_yerr = np.hypot(e600, e750) if np.isfinite(e600) and np.isfinite(e750) else None
    selected_mark(axes[0, 1], source_phot('m475') - source_phot('m600'), source_phot('m600') - source_phot('m750'), ccd_xerr, ccd_yerr)
    axes[0, 1].set(xlabel='m475 - m600', ylabel='m600 - m750', title='Color-color diagram')
    if axes[0, 1].get_legend_handles_labels()[0]: axes[0, 1].legend(fontsize=8)
    mag = phot(population_catalog, diagnostic_band); err = phot(population_catalog, diagnostic_band, 'mag_err')
    snr = np.divide(1.0857, err, out=np.full(len(err), np.nan), where=np.isfinite(err) & (err > 0))
    valid_sn = np.isfinite(mag) & np.isfinite(snr)
    axes[0, 2].scatter(mag[valid_sn], snr[valid_sn], s=6, alpha=0.16, color='black', rasterized=True)
    if valid_sn.any():
        edges = np.arange(np.floor(np.nanmin(mag[valid_sn])), np.ceil(np.nanmax(mag[valid_sn])) + 0.25, 0.25)
        centers, medians = [], []
        for left, right in zip(edges[:-1], edges[1:]):
            sel = valid_sn & (mag >= left) & (mag < right)
            if np.count_nonzero(sel) >= 10: centers.append((left + right) / 2); medians.append(np.nanmedian(snr[sel]))
        axes[0, 2].plot(centers, medians, 'o-', color='tab:blue', label='bin median (N >= 10)')
        sn5 = _sn5_limit(centers, medians)
        if np.isfinite(sn5): axes[0, 2].axvline(sn5, color='tab:orange', linestyle='-.', label=f'S/N=5 mag: {sn5:.2f}')
    selected_error = source_phot(diagnostic_band, 'mag_err')
    selected_mark(axes[0, 2], source_phot(diagnostic_band), 1.0857 / selected_error if np.isfinite(selected_error) else np.nan)
    axes[0, 2].axhline(5, color='black', linestyle='--', linewidth=1, label='S/N = 5')
    axes[0, 2].set(xlabel=f'{diagnostic_band} magnitude', ylabel='S/N', title=f'{diagnostic_band} S/N-Mag'); axes[0, 2].legend(fontsize=8)
    hist_values = mag[np.isfinite(mag)]
    if hist_values.size:
        selected_mag = source_phot(diagnostic_band)
        hist_lo, hist_hi = np.nanpercentile(hist_values, [0.5, 99.5])
        if np.isfinite(selected_mag): hist_lo, hist_hi = min(hist_lo, selected_mag), max(hist_hi, selected_mag)
        hist_edges = np.arange(np.floor(hist_lo), np.ceil(hist_hi) + 0.25, 0.25)
        axes[1, 0].hist(hist_values, bins=hist_edges, histtype='step', color='black')
        if np.isfinite(selected_mag):
            bin_index = min(np.searchsorted(hist_edges, selected_mag, side='right') - 1, len(hist_edges) - 2)
            if 0 <= bin_index < len(hist_edges) - 1:
                axes[1, 0].axvspan(hist_edges[bin_index], hist_edges[bin_index + 1], color='red', alpha=0.25)
                bin_count = int(np.count_nonzero((hist_values >= hist_edges[bin_index]) & (hist_values < hist_edges[bin_index + 1])))
                axes[1, 0].axvline(selected_mag, color='red', label=f'source bin: {hist_edges[bin_index]:.2f}-{hist_edges[bin_index + 1]:.2f} (N={bin_count})')
                axes[1, 0].legend()
    axes[1, 0].set(xlabel=f'{diagnostic_band} magnitude', ylabel='Number of objects', title='Magnitude histogram')
    sed_mags = np.asarray(mags, dtype=float)
    sed_mag_errors = np.asarray(magerrs, dtype=float)
    sed_valid = np.isfinite(sed_mags) & np.isfinite(sed_mag_errors) & (sed_mag_errors > 0)
    # Magnitudes are AB; converting +/- magnitude errors gives asymmetric flux errors.
    sed_flux_ujy = 3631.0e6 * 10.0 ** (-0.4 * sed_mags)
    sed_flux_faint_ujy = 3631.0e6 * 10.0 ** (-0.4 * (sed_mags + sed_mag_errors))
    sed_flux_bright_ujy = 3631.0e6 * 10.0 ** (-0.4 * (sed_mags - sed_mag_errors))
    sed_flux_errors_ujy = np.vstack((sed_flux_ujy - sed_flux_faint_ujy, sed_flux_bright_ujy - sed_flux_ujy))
    sed_valid &= np.isfinite(sed_flux_ujy) & (sed_flux_ujy > 0)
    axes[1, 1].set(xlabel='Central wavelength (nm)', ylabel=r'$f_\nu$ ($\mu$Jy)', title=f'7DT {N_ACTIVE_BANDS}-band flux-density SED')
    axes[1, 1].set_xlim(375, 900)
    if np.any(sed_valid):
        axes[1, 1].errorbar(np.asarray(wavelengths)[sed_valid], sed_flux_ujy[sed_valid], yerr=sed_flux_errors_ujy[:, sed_valid], fmt='o', color='tab:purple')
        axes[1, 1].set_yscale('log')
    else:
        axes[1, 1].text(0.5, 0.5, 'No valid 7DT photometry', ha='center', va='center', transform=axes[1, 1].transAxes)
    # Summary values are collected once and used by both screen and PDF output.
    summary['COSMOS'] = [COSMOS]
    if unified:
        summary['field_master_id'] = [int(original_row['selected_field_master_id'])]
        summary['n_success_bands_across_fields'] = [int(original_row['n_success_bands'])]
        summary['n_available_bands_across_fields'] = [int(original_row['n_available_bands'])]
        if 'finalcat_RM' in original_row.colnames:
            summary['finalcat_RM'] = [original_row['finalcat_RM']]
    _style_image_page(report_fig, image_ax, band_axes, axes)
    try:
        table_fig = _table_page(summary, band_table, finalcat_table, sdss_table, source_id, COSMOS)
    except Exception:
        plt.close(report_fig)
        raise
    result = SourceInfoResult({'master_id': selected_id, 'row': row,
        'image': (image_fig, image_ax), 'summary_table': summary, 'cutout_figure': band_fig,
        'band_table': band_table, 'finalcat_table': finalcat_table, 'sdss_table': sdss_table,
        'diagnostic_figure': pop_fig, 'report_figures': (report_fig, table_fig), 'pdf_path': None})
    if save_pdf is not False and save_pdf is not None:
        if save_pdf is True:
            directory = Path(output_dir) if output_dir is not None else project_root / f'output/COSMOS{COSMOS}/singlephot/analysis/source_info'
            local = f"_master{int(original_row['selected_field_master_id'])}" if unified else ''
            pdf_path = directory / f'cosmos{COSMOS}_source{selected_id}{local}_sourceinfo.pdf'
        elif isinstance(save_pdf, (str, Path)):
            pdf_path = Path(save_pdf)
        else:
            raise TypeError('save_pdf must be False, True, or a PDF path.')
        pdf_path = pdf_path.expanduser().resolve()
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with PdfPages(pdf_path, metadata={'Title': f'COSMOS{COSMOS} source {selected_id}'}) as pdf:
                for figure in result['report_figures']:
                    # Do not use bbox_inches='tight': physical A4 page dimensions must stay exact.
                    pdf.savefig(figure, dpi=200)
        finally:
            if not show:
                for figure in result['report_figures']:
                    plt.close(figure)
        result['pdf_path'] = pdf_path
    if show:
        # Explicit display also works when Agg is used to validate notebook code.
        from IPython.display import display
        for figure in result['report_figures']:
            display(figure)
            plt.close(figure)
    else:
        for figure in result['report_figures']:
            plt.close(figure)
    return result
def plot_source_image(cat, image=None, pixel_range=None, source_id=None, ax=None):
    '''Plot a detection-image cutout centered on a selected catalog source.'''
    if image is None:
        raise ValueError('A detection image is required.')
    if len(cat) == 0:
        raise ValueError('The input catalog is empty.')
    required = {'master_id', 'master_x', 'master_y'}
    missing = required.difference(cat.colnames)
    if missing:
        raise KeyError(f'Catalog is missing required columns: {sorted(missing)}')
    if pixel_range is not None and (not np.isfinite(pixel_range) or pixel_range <= 0):
        raise ValueError('pixel_range must be None or a positive number.')
    master_ids = np.asarray(cat['master_id'], dtype=np.int64)
    if source_id is None:
        selected_index = int(np.random.default_rng().integers(len(cat)))
    else:
        matches = np.flatnonzero(master_ids == int(source_id))
        if len(matches) == 0:
            raise KeyError(f'master_id={source_id} is not present in the catalog.')
        selected_index = int(matches[0])
    x_all = np.ma.filled(cat['master_x'], np.nan).astype(float)
    y_all = np.ma.filled(cat['master_y'], np.nan).astype(float)
    selected_x, selected_y = float(x_all[selected_index]), float(y_all[selected_index])
    if not np.isfinite(selected_x) or not np.isfinite(selected_y):
        raise ValueError('The selected source has a non-finite pixel position.')
    if pixel_range is None:
        xmin, xmax, ymin, ymax = 0, image.shape[1], 0, image.shape[0]
    else:
        center_x, center_y = int(round(selected_x - 1)), int(round(selected_y - 1))
        half = int(round(pixel_range))
        xmin, xmax = max(0, center_x - half), min(image.shape[1], center_x + half + 1)
        ymin, ymax = max(0, center_y - half), min(image.shape[0], center_y + half + 1)
    finite = np.isfinite(x_all) & np.isfinite(y_all)
    in_view = finite & (x_all - 1 >= xmin) & (x_all - 1 < xmax) & (y_all - 1 >= ymin) & (y_all - 1 < ymax)
    crop = image[ymin:ymax, xmin:xmax]
    fig = ax.figure
    ax.imshow(crop, origin='lower', cmap='gray', norm=ImageNormalize(image, interval=ZScaleInterval()), extent=[xmin + 0.5, xmax + 0.5, ymin + 0.5, ymax + 0.5])
    class_colors = {'matched': 'purple', '7dt_only': 'tab:red', 'finalcat_only': 'tab:blue'}
    classes = np.asarray(cat['master_class']).astype(str) if 'master_class' in cat.colnames else np.full(len(cat), 'source')
    for source_class in np.unique(classes[in_view]):
        selected = in_view & (classes == source_class)
        ax.scatter(x_all[selected], y_all[selected], s=24, facecolors='none', edgecolors=class_colors.get(source_class, 'tab:red'), linewidths=0.8, label=f'{source_class} ({int(selected.sum())})')
    ax.scatter(selected_x, selected_y, marker='+', s=180, color='cyan', linewidths=1.5, label=f'master_id={master_ids[selected_index]}')
    ax.set(xlabel='X pixel', ylabel='Y pixel')
    title_class = classes[selected_index]
    view_label = 'whole image' if pixel_range is None else f'pixel_range={pixel_range}'
    ax.set_title(f'{title_class}: master_id={master_ids[selected_index]} ({view_label})')
    ax.legend(loc='upper right')
    return fig, ax, int(master_ids[selected_index])


def _sn5_limit(centers, medians, level=5.0):
    # Median S/N falls with magnitude; linearly interpolate the bright-to-faint
    # crossing of the S/N = level line. Returns nan when there is no crossing.
    c = np.asarray(centers, dtype=float)
    m = np.asarray(medians, dtype=float)
    ok = np.isfinite(c) & np.isfinite(m)
    c, m = c[ok], m[ok]
    if c.size < 2:
        return np.nan
    order = np.argsort(c)
    c, m = c[order], m[order]
    below = np.where(m < level)[0]
    if below.size == 0 or np.all(m < level):
        return np.nan
    i = int(below[0])
    if i == 0 or not (m[i - 1] >= level > m[i]):
        return np.nan
    c0, m0, c1, m1 = c[i - 1], m[i - 1], c[i], m[i]
    return float(c0 + (level - m0) * (c1 - c0) / (m1 - m0))
