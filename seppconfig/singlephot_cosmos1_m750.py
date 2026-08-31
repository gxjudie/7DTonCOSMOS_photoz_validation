
from sourcextractor.config import *

image_group = load_fits_image(
    '/Users/mac_seo/Research_local/7DTonCOSMOS/COSMOS_data/COSMOS1/calib_7DT_COSMOS_1_20240511_060100_m750_8760.com.fits',
    # weight='/Users/mac_seo/Research_local/7DTonCOSMOS/COSMOS_data/COSMOS1/calib_7DT_COSMOS_1_20240511_060100_m750_8760.weight.fits',
    # weight_type='weight'
)

meas_group = MeasurementGroup(image_group)

aper = []

for img in meas_group:
    aper.extend(
        add_aperture_photometry(
            img,
            5.940594059405
        )
    )

add_output_column("aper3", aper)
