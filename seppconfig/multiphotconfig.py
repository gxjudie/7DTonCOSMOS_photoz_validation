from sourcextractor.config import *
from pathlib import Path

# ---------------------------------------------------------
# 1. Initial Setting
# ---------------------------------------------------------
COSMOS = 1
DATA_ROOT = Path("/Users/mac_seo/Research_local/7DTonCOSMOS/COSMOS_data")
BANDS_7DT = sorted({f'm{i}' for i in range(400, 900, 25)})
pixel_3arcsecaper = 5.94059406

def get_image_path(cosmos, band, weight=False):
    image_folder = DATA_ROOT / f"COSMOS{cosmos}"
    suffix = "weight.fits" if weight else "com.fits"
    paths = list(image_folder.glob(f"*{band}*.{suffix}"))

    if len(paths) == 0:
        return None
    return str(paths[0])

image_paths = []
for i, band in enumerate(BANDS_7DT):
    path = get_image_path(COSMOS, band)
    if path:
        image_paths.append(path)

# ---------------------------------------------------------
# 2. Photometry
# ---------------------------------------------------------

image_group = load_fits_images(image_paths)
meas_group = MeasurementGroup(image_group)
 
aper = [] 
 
for img in meas_group : 
    aper.extend(add_aperture_photometry(img, pixel_3arcsecaper))

add_output_column("aper3", aper)

