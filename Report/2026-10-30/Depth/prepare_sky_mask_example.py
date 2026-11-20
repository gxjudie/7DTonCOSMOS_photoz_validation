"""Render the presentation's COSMOS1 m650 masking example from frozen inputs.

Only the display changes. The raw centroids, common valid mask and archived
flag-cut size proxy retain the depth analysis's original definitions. This
script creates one figure and never rewrites scientific data products.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
FIELD, BAND = 1, 'm650'
BASE = ROOT / 'output/experiments/aperture_flag_ablation/baseline/full_products'


def main():
    """Compare identical full-resolution central-cell cutouts before/after masking."""
    raw = Table.read(ROOT / f'output/COSMOS{FIELD}/singlephot/cat_cosmos{FIELD}_{BAND}.fits', hdu=1)
    xy = np.column_stack([
        np.ma.filled(raw['pixel_centroid_x'].astype(float), np.nan) - 1,
        np.ma.filled(raw['pixel_centroid_y'].astype(float), np.nan) - 1,
    ])
    # Mask all finite raw centers, including detections rejected from the
    # final photometry sample; the FWHM scale uses the separate baseline table.
    xy = xy[np.isfinite(xy).all(axis=1)]
    sizes = Table.read(BASE / f'COSMOS{FIELD}/singlephot/diagnostics/cosmos{FIELD}_7dt_singleband_band_uncertainty.fits')
    labels = np.char.strip(np.asarray(sizes['band']).astype(str))
    fwhm = float(sizes['fwhm_pix'][labels == BAND][0])
    radius = 5 * fwhm / 2.355
    image_path = next((ROOT / f'COSMOS_data/COSMOS{FIELD}').glob(f'*{BAND}*.com.fits'))
    with fits.open(image_path, memmap=True) as hdus:
        ny, nx = hdus[0].data.shape
        # Show one central 10x10-grid cell without the old stride-32 sampling,
        # which made the circular source masks too small to see on a slide.
        row, column = 4, 4
        x0, x1 = nx * column // 10, nx * (column + 1) // 10
        y0, y1 = ny * row // 10, ny * (row + 1) // 10
        cutout = np.asarray(hdus[0].data[y0:y1, x0:x1], dtype=float)
        pixel_scale = float(proj_plane_pixel_scales(WCS(hdus[0].header))[0] * 3600)
    with fits.open(ROOT / f'COSMOS_data/COSMOS{FIELD}/COSMOS{FIELD}_validmask.fits', memmap=True) as hdus:
        valid = np.asarray(hdus[0].data[y0:y1, x0:x1], dtype=bool)
    yy, xx = np.indices(cutout.shape)
    distance, _ = cKDTree(xy).query(np.column_stack([(xx+x0).ravel(), (yy+y0).ravel()]))
    source_mask = distance.reshape(cutout.shape) <= radius
    masked = np.ma.masked_where(source_mask | ~valid | ~np.isfinite(cutout), cutout)
    low, high = ZScaleInterval(contrast=.2).get_limits(cutout[np.isfinite(cutout)])
    norm = ImageNormalize(vmin=low, vmax=high, stretch=AsinhStretch(.3))
    cmap = plt.get_cmap('gray').copy()
    cmap.set_bad('#D57C44')
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), sharex=True, sharey=True)
    for ax, data, title in zip(axes, [cutout, masked], ['Original science image', 'Source-masked sky']):
        ax.imshow(data, origin='lower', cmap=cmap, norm=norm,
                  extent=(x0-.5, x1-.5, y0-.5, y1-.5), interpolation='nearest')
        ax.set_title(title, fontsize=15)
        ax.set_xlabel('Original image X [pixel]')
    axes[0].set_ylabel('Original image Y [pixel]')
    fig.suptitle(f'COSMOS1 m650: circular source masks, r = {radius:.2f} pix ({radius*pixel_scale:.2f} arcsec)', fontsize=16)
    fig.text(.07, .025, 'Orange: excluded pixels. Same cutout and intensity scale. All finite raw centroids define mask centers.', fontsize=10)
    fig.tight_layout(rect=(0, .13, 1, .95))
    target = OUT / 'cosmos1_7dt_depth_m650_sky_mask_example.png'
    fig.savefig(target, dpi=150)
    plt.close(fig)
    print(target)


if __name__ == '__main__':
    main()
