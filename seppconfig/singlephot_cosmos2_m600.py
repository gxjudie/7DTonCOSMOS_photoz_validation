
from sourcextractor.config import *

image_group = load_fits_image('/Users/mac_seo/Research_local/7DTonCOSMOS/COSMOS_data/COSMOS2/calib_7DT_COSMOS_2_20240527_234344_m600_4680.com.fits')
meas_group = MeasurementGroup(image_group)
aper = []
for img in meas_group:
    aper.extend(add_aperture_photometry(img, 5.940594059405))
add_output_column("aper3", aper)
