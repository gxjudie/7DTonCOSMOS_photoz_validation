
from sourcextractor.config import *

image_group = load_fits_image('/Users/mac_seo/Research_local/7DTonCOSMOS/COSMOS_data/COSMOS3/calib_7DT_COSMOS_3_20240526_000116_m850_2760.com.fits')
meas_group = MeasurementGroup(image_group)
aper = []
for img in meas_group:
    aper.extend(add_aperture_photometry(img, 5.940594059405))
add_output_column("aper3", aper)
