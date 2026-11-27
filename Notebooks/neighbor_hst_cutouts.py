"""Cache native-resolution HST/ACS COSMOS display cutouts for source inspection.

This module queries images only. Source associations must already exist in the
7DT master; no catalog matching or source reassignment is performed here.
"""
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import json
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

SIA_URL = 'https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas'


def fetch_hst_cutout(cache_dir, source_id, ra, dec, width_arcsec, *, cache_prefix="m600"):
    """Return (FITS path, coverage fraction), or None for absent ACS coverage.

    ``source_id`` is a cache identifier, not a new source association.
    Use ``cache_prefix="finalcat"`` when the identifier is a reference row key.
    Reproject the overlapping original ACS F814W science tiles to a common TAN
    grid at approximately 0.03 arcsec/pixel. The unweighted overlap average is
    for morphology display, not aperture photometry. All input URLs are saved.
    Network failures propagate instead of being mistaken for missing coverage.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f'{cache_prefix}_{source_id}_hst_acs_f814w.fits'
    meta_path = path.with_suffix('.json')
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        assert np.allclose([meta['ra'], meta['dec'], meta['width_arcsec']],
                           [ra, dec, width_arcsec], rtol=0, atol=1.e-8)
        if path.exists():
            return path, meta['coverage_fraction']
        if meta.get('no_coverage'):
            return None
    response = requests.get(SIA_URL, params={
        'mission': 'COSMOS', 'hdr_location': '\\COSMOSDataPath\\',
        'collection_desc': 'Cosmic Evolution Survey with HST (COSMOS)',
        'SIAP_ACTIVE': '1', 'POS': f'{ra},{dec}',
        'SIZE': str((width_arcsec + 8) / 3600), 'FORMAT': 'image/fits'}, timeout=45)
    response.raise_for_status()
    table = parse_single_table(BytesIO(response.content)).to_table()
    # Restrict to the original F814W ACS science tiles; no weights, previews,
    # synthetic images, or mixed-filter mosaics are substituted.
    tiles = sorted({str(row['col35']) for row in table
                    if str(row['col29']).lower() == 'science'
                    and str(row['col35']).startswith('images/acs_2.0/I/')})
    meta = {'source_id': int(source_id), 'cache_prefix': cache_prefix, 'ra': float(ra), 'dec': float(dec),
            'width_arcsec': float(width_arcsec), 'sia_url': response.url}
    if not tiles:
        meta.update(no_coverage=True, coverage_fraction=0.)
        meta_path.write_text(json.dumps(meta, indent=2)+'\n')
        return None
    npix = int(np.ceil(width_arcsec / .03))
    output_wcs = WCS(naxis=2)
    output_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    output_wcs.wcs.crval = [ra, dec]
    output_wcs.wcs.crpix = [(npix+1)/2, (npix+1)/2]
    output_wcs.wcs.cdelt = [-width_arcsec/3600/npix, width_arcsec/3600/npix]
    yy, xx = np.indices((npix, npix), dtype=float)
    world_ra, world_dec = output_wcs.pixel_to_world_values(xx, yy)
    total = np.zeros((npix, npix))
    count = np.zeros((npix, npix), dtype=np.uint8)
    urls = []
    for tile in tiles:
        cutout = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/Subimage/nph-subimage',
                              params={'origfile': '/irsadata/COSMOS/'+tile,
                                      'ra': str(ra), 'dec': str(dec),
                                      'xsize': str((width_arcsec+8)/3600)}, timeout=45)
        cutout.raise_for_status()
        urls.append(cutout.url)
        with fits.open(BytesIO(cutout.content), memmap=False) as hdus:
            hdu = next(h for h in hdus if h.data is not None and h.data.ndim == 2)
            data = np.asarray(hdu.data, dtype=float).copy()
            image_wcs = WCS(hdu.header)
        # The rotated ACS tile bounding boxes have zero-valued padding.
        data[data == 0] = np.nan
        px, py = image_wcs.world_to_pixel_values(world_ra, world_dec)
        projected = map_coordinates(data, [py, px], order=1, mode='constant',
                                    cval=np.nan, prefilter=False)
        good = np.isfinite(projected)
        total[good] += projected[good]
        count[good] += 1
    mosaic = np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)
    coverage = float(np.mean(count > 0))
    header = output_wcs.to_header()
    header['SOURCEID'] = int(source_id)
    header['COVERAGE'] = coverage
    fits.PrimaryHDU(mosaic.astype(np.float32), header=header).writeto(path, overwrite=True)
    meta.update(input_urls=urls, coverage_fraction=coverage, no_coverage=False,
                display_only=True, projection='north-up TAN; bilinear; mean overlap')
    meta_path.write_text(json.dumps(meta, indent=2)+'\n')
    return path, coverage


def select_hst_covered(rows, fetch_one, needed=9, batch_size=4):
    """Inspect a fixed candidate order, choosing only >=99% covered cutouts.

    Batched downloads preserve candidate ordering despite parallel completion.
    Errors are raised; only documented absence/partial coverage skips a source.
    """
    selected, images = [], []
    start = 0
    while start < len(rows):
        take = min(batch_size, needed-len(selected))
        batch = rows[start:start+take]
        start += len(batch)
        with ThreadPoolExecutor(max_workers=batch_size) as executor:
            results = list(executor.map(fetch_one, batch))
        for row, result in zip(batch, results):
            if result is not None and result[1] >= .99:
                selected.append(row)
                images.append(result[0])
                if len(selected) == needed:
                    return np.asarray(selected, dtype=int), images
    raise ValueError(f'Only {len(selected)} eligible sources have >=99% ACS coverage.')
