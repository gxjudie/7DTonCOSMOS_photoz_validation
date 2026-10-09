# 7DT Single-band Matching in COSMOS1

## Slide 1 — Input and workflow (about 45 seconds)

Hello. Today I will show how I matched sources in the 7DT COSMOS1 data. This catalog is one step toward testing redshift estimates. I started with 20 separate band catalogs. I found and measured sources in each band. After quality cuts, I used 78,882 detections from the shared area. I first linked the 7DT detections across bands. Then I compared the new 7DT catalog with FINALCAT.

## Slide 2 — 7DT cross-band matching (about 65 seconds)

The images share the same pixel grid. I processed the bands in wavelength order, starting with m400. A source position is less certain when the image is less sharp or the signal is weaker. So I gave each band its own position error. For each new source, I compared its distance from the current master with their combined error. I kept pairs within three times that error. When a source joined a master, I updated the master position. Better positions received more weight. This made one position for each linked object.

## Slide 3 — Internal check (about 65 seconds)

Next, I checked the final 7DT masters. For each member, I measured its distance from the final master position. I divided this distance by the expected position error. Most members are close to their final master. In total, 99.85 percent are within three times the error. Only 117 of 78,882 members are outside this range. The final catalog has 19,489 7DT masters. This check shows good position agreement within the catalog. It does not prove that every match is correct.

## Slide 4 — FINALCAT radius (about 75 seconds)

I then matched the 7DT masters to FINALCAT using sky positions. I tested three maximum distances. To estimate chance matches, I shifted the FINALCAT positions and ran the match again. I used eight shifted samples. At 0.5 arcseconds, there were 2,682 real matches. At 0.7 arcseconds, there were 2,977. The share of chance matches was about 3.6 percent. At 1.0 arcseconds, the number of real matches rose to 3,212. But the share of chance matches rose faster, to about 6.6 percent. I therefore chose 0.7 arcseconds. The median distance of the selected pairs is 0.205 arcseconds.

## Slide 5 — Final catalog (about 50 seconds)

The final result has 2,977 objects seen in both 7DT and FINALCAT. It also has 16,512 objects found only in 7DT, and 55,770 found only in FINALCAT. Together, they form a union catalog of 75,259 objects. The main points are simple. The 7DT matches use the position error of each band. The internal position check is strong. The FINALCAT radius has a chance-match test behind it. This catalog is ready for the next photo-z analysis. Thank you.
