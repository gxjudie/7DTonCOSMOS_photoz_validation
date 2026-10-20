"""Reuse complete SDSS candidate caches without changing matching selection.

A candidate query is replaced only when its entire sky-search circle fits
inside at least one cached DR17 primary-object SQL rectangle. Return every
nearby candidate, including ambiguous ones; callers retain their normal
one-to-one association and conflict-resolution rules.
"""
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table


def complete_sdss_candidates(row_indices, reference_coordinates, cache_paths,
                             column_names, radius_arcsec, payloads_equal):
    """Return all cached candidates, or None when completeness is unproven."""
    targets=reference_coordinates[row_indices]
    radius=float(radius_arcsec)*u.arcsec
    dec=np.asarray(targets.dec.deg)
    ra=np.asarray(targets.ra.deg)
    dec_margin=radius.to_value(u.deg)
    ra_margin=dec_margin/np.cos(np.deg2rad(np.abs(dec)+dec_margin))
    covered=np.zeros(len(row_indices),bool)
    catalogs=[]
    for cache_path in cache_paths:
        if not Path(cache_path).exists():
            continue
        table=Table.read(cache_path)
        if (int(table.meta.get('SDSS_DR',-1))!=17 or not set(column_names).issubset(table.colnames)
                or not np.all(np.asarray(table['mode'])==1)):
            continue
        meta=table.meta
        covered|=((ra-ra_margin>=float(meta.get('RA_MIN',np.inf)))
                  & (ra+ra_margin<=float(meta.get('RA_MAX',-np.inf)))
                  & (dec-dec_margin>=float(meta.get('DEC_MIN',np.inf)))
                  & (dec+dec_margin<=float(meta.get('DEC_MAX',-np.inf))))
        catalogs.append(table)
    if not covered.all():
        return None
    payloads={}
    for table in catalogs:
        coordinates=SkyCoord(np.asarray(table['ra'],float)*u.deg,
                             np.asarray(table['dec'],float)*u.deg)
        _,selected,_,_=search_around_sky(targets,coordinates,radius)
        for index in np.unique(selected):
            row=table[int(index)]
            payload={name:row[name] for name in column_names}
            objid=int(row['objid'])
            if objid in payloads and not payloads_equal(payloads[objid],payload):
                return None
            payloads[objid]=payload
    print(f'SDSS candidate retrieval: {len(row_indices)} complete search circles reused from local DR17 caches.')
    return payloads
