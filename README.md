# Validation of 7DT Photometric Redshift Performance in the COSMOS Field

This repository contains an ongoing research project to validate the photometric redshift (photo-z) performance of the 7-Dimensional Telescope (7DT) using real COSMOS-field observations.

The main goal is to test whether the 7DT medium-band filter system can produce reliable photo-z estimates from observed galaxy photometry, and to assess whether the resulting catalog may be useful for future weak-lensing applications.

## Research Motivation

Photometric redshifts are essential for large imaging surveys because spectroscopic redshifts are too expensive to obtain for every faint extragalactic source. However, photo-z estimates can suffer from scatter, systematic bias, and catastrophic outliers. These errors are especially important for weak-lensing studies, where reliable source redshift distributions are required.

The 7DT medium-band filter system provides denser wavelength sampling than conventional broad-band surveys. This should help trace galaxy spectral energy distributions (SEDs), including features such as the 4000 Angstrom break and strong emission lines. Previous work forecasted strong photo-z performance for 7DT using mock catalogs. This project aims to test that forecast with real 7DT COSMOS observations.

## Objectives

1. Construct a real multi-band photometric catalog from 7DT observations of the COSMOS field.
2. Validate 7DT photo-z performance by comparing estimated photo-z values with external reference redshifts.
3. Assess whether the resulting catalog is suitable for future weak-lensing applications.

## Data

The project uses calibrated 7DT FITS images of the COSMOS 1-4 fields and corresponding filter response curves.

The 7DT filter set used in this study consists of 20 medium-band filters:

`m400`, `m425`, `m450`, `m475`, `m500`, `m525`, `m550`, `m575`, `m600`, `m625`, `m650`, `m675`, `m700`, `m725`, `m750`, `m775`, `m800`, `m825`, `m850`, `m875`

The number in each filter name indicates the approximate central wavelength in nanometers.

Large FITS data files and generated catalogs are not tracked in git. They are expected to exist locally under `COSMOS_data/` and `output/`.

The nominal filter set contains 20 bands, but each field's active bands are determined from the available science/weight image pairs. COSMOS4 uses binary coverage weights created by `0.1_make_COSMOS4_boolean_weights.ipynb`; these encode finite, nonzero science pixels and are not inverse-variance maps.

## Method Summary

The intended analysis workflow is:

1. Inspect calibrated 7DT images and WCS consistency.
2. Build valid masks for regions covered by all usable bands.
3. Generate single-band and/or multi-band detection catalogs using SourceXtractor++.
4. Match single-band catalogs, associate valid-footprint FINALCAT sources, and enrich each field's master catalog with SDSS DR17 Petrosian photometry.
5. Combine FINALCAT-matched sources across COSMOS1–4, inspect band coverage and depth, and prepare multi-band SEDs for fitting.
6. Perform SED fitting using CIGALE with 7DT transmission curves.
7. Estimate photometric redshifts from the 7DT-based photometry.
8. Compare photo-z estimates with spectroscopic or high-quality external reference redshifts.
9. Quantify performance using catastrophic failure fraction, sigma_NMAD, and bias.

## Photo-z Validation Metrics

For reference redshift `z_ref` and photometric redshift `z_phot`, the residual is defined as:

```text
Delta z = z_ref - z_phot
```

The catastrophic failure fraction is defined by:

```text
eta = fraction(|z_ref - z_phot| / (1 + z_ref) > 0.15)
```

The scatter is measured using normalized median absolute deviation:

```text
sigma_NMAD = 1.48 * median(|Delta z - median(Delta z)| / (1 + z_ref))
```

The systematic bias is defined as:

```text
bias = mean((z_phot - z_ref) / (1 + z_ref))
```

## Repository Structure

```text
.
├── Notebooks/
│   ├── 0_preReaserch.ipynb
│   ├── 0.1_make_COSMOS4_boolean_weights.ipynb
│   ├── 1_multibandphot.ipynb
│   ├── 2_singphotandmatching.ipynb
│   ├── 2.1_singlephotcatalog_visualization.ipynb
│   ├── 2.2_7DTdepthAnalysis.ipynb
│   ├── 2.3_numofmatchedbands_digonostics.ipynb
│   ├── sourceinfo_report.py
│   └── 3_photozvalidation.ipynb
├── seppconfig/
│   ├── multiphotconfig.config
│   ├── multiphotconfig.py
│   ├── multiphotconfig_cosmos{2,3,4}.config / .py
│   └── singlephot_cosmos{1,2,3,4}_m*.config / .py
├── sed_fitting/
│   ├── pcigale.ini
│   ├── pcigale.ini.spec
│   └── transmission_curve/
├── 7DTonCOSMOS_tex/
│   └── Proposal/
├── Report/             # tracked input CSVs, report text, and build scripts
│                       # generated figures/PDFs/slides stay local
├── COSMOS_data/        # local data, ignored by git
├── output/             # generated products, ignored by git
└── cigale/             # local/vendor CIGALE copy, ignored by git
```

## Notebook Roles

- `0_preReaserch.ipynb`: preliminary checks, WCS consistency tests, valid-mask construction, S/N checks, and inspection of the external COSMOS catalog.
- `0.1_make_COSMOS4_boolean_weights.ipynb`: create COSMOS4 binary coverage-weight images while retaining image WCS information.
- `1_multibandphot.ipynb`: multi-band photometry experiments, chi/chi-mean detection-image construction, SourceXtractor++ execution, and SED plot checks.
- `2_singphotandmatching.ipynb`: production single-band photometry, cross-band matching, FINALCAT union, SDSS integration, completeness curves, and positional-quality products for a selected field.
- `2.1_singlephotcatalog_visualization.ipynb`: read-only visual analysis of the authoritative catalog products.
- `2.2_7DTdepthAnalysis.ipynb`: source inspection, magnitude/S/N and size diagnostics, blank-aperture depth maps, combined COSMOS1–4 depth maps, and external-catalog recovery analysis.
- `2.3_numofmatchedbands_digonostics.ipynb`: construct and reload a FINALCAT-centered COSMOS1–4 catalog, reconcile SDSS counterparts, and analyze successful band counts and FINALCAT RM distributions.
- `sourceinfo_report.py`: shared source viewer for notebooks 2.2 and 2.3, with two portrait A4 pages containing image cutouts, diagnostics, SED, and source tables; supports PDF export.
- `3_photozvalidation.ipynb`: placeholder for final photo-z validation; currently empty.

## Running the Catalog and Diagnostics Workflow

Use a Python/Jupyter environment with NumPy, Astropy, SciPy, Matplotlib, and tqdm. SourceXtractor++ must be available as `sourcextractor++` when generating missing single-band catalogs. SDSS downloads and conflict rechecks require network access. CIGALE is needed for the later SED-fitting stage.

1. Place science images (`*.com.fits`), weight images (`*.weight.fits`), and `COSMOS{field}_validmask.fits` under `COSMOS_data/COSMOS{field}/`; place the reference catalog at `COSMOS_data/COSMOS_finalcat.fits`. Create COSMOS4 binary weights if needed.
2. Run `2_singphotandmatching.ipynb` with `COSMOS` set to each field from 1 through 4. Review `OVERWRITE_PHOTOMETRY`, `OVERWRITE_MATCHING`, and `OVERWRITE_SDSS`: the current defaults reuse existing raw photometry but rebuild matching and SDSS products. The notebook generates field/band SourceXtractor++ configuration files when photometry runs.
3. Use `2.1_singlephotcatalog_visualization.ipynb` for visual inspection of saved production products. Run the per-field depth sections in `2.2_7DTdepthAnalysis.ipynb` before its combined-field depth maps, which consume the saved depth grids.
4. Run `2.3_numofmatchedbands_digonostics.ipynb` after all four field master catalogs exist. It uses their existing FINALCAT links, writes the unified catalog, and reloads it for subsequent diagnostics.

Launch notebooks from the project root or `Notebooks/` so their project-root detection resolves local data correctly. The multi-band detection experiments remain in notebook 1; the production single-band catalog workflow is in notebook 2.

## Catalog Definitions and Products

Per-field products are stored in `output/COSMOS{field}/singlephot/`:

| Product | Location relative to the field's `singlephot/` directory |
| --- | --- |
| Raw single-band catalogs | `cat_cosmos{field}_{band}.fits` |
| Authoritative SDSS-enriched master | `cosmos{field}_master_catalog.fits` |
| Pre-SDSS 7DT–FINALCAT union | `diagnostics/cosmos{field}_7dt_finalcat_matching_union_catalog.fits` |
| Complete available-band subset | `diagnostics/cosmos{field}_master_population_fullband_catalog.fits` |
| Matching, positional-quality, and population diagnostics | `diagnostics/` |
| Recovery curves, completeness depths, and blank-aperture depth products | `analysis/` |

The field master contains `matched`, `7dt_only`, and `finalcat_only` populations. FINALCAT associations use a 0.7-arcsec radius; independent SDSS DR17 primary-object matches use a deterministic one-to-one assignment within 1.5 arcsec and do not add master rows. SDSS magnitudes and errors are stored as `sdss_petroMag_{u,g,r,i,z}` and `sdss_petroMagErr_{u,g,r,i,z}`. FINALCAT-provided SDSS information remains under `finalcat_SDSS_*`.

The combined catalog is written to:

```text
output/singlephot/master_catalog_MatchedWithFianlcat_OverallCosmosField.fits
```

The spelling `Fianlcat` is part of the current filename. This catalog contains one row per `finalcat_idx` with a 7DT match in at least one field. The original FINALCAT ID string is not assumed to be unique. Shared reference information is retained as `finalcat_*`, field measurements as `7dt_cosmos{field}_{original_name}`, and reconciled SDSS information as common `sdss_*` columns. Global `master_id` values are separate from the original field-local IDs. Conflicting SDSS matches are rechecked; unresolved ambiguous counterparts remain masked, with decisions recorded in CSV logs.

The unified FITS contains `CATALOG`, `COLUMN_MASKS`, `FIELD_METADATA`, and `COLUMN_MAP` extensions. Use `read_overall_catalog` from notebook 2.3 to restore all column masks, including Boolean and string masks.

For notebook 2.3, a band succeeds when its source ID is unmasked and nonnegative, its magnitude is finite with `0 < mag < 50`, and its magnitude error is finite and positive. No additional S/N cut is applied. Success is combined across matched fields, counting each distinct band once. `n_success_bands` therefore ranges from 0 to 20; `all_20_success` requires all 20 nominal filters, while `all_available_success` requires all filters available across that source's matched fields. This differs from the per-field complete available-band subset, which uses production detection membership.

Band-count histograms, RM distributions, the ≥15-band comparison, SDSS reconciliation logs, and `summary.json` are saved under `output/singlephot/analysis/numofmatchedbands/`. `finalcat_RM` is the FINALCAT Subaru r+ total AB magnitude. Source inspection uses a representative field chosen by successful-band count, then median magnitude error, FINALCAT separation, and field number; cutouts and photometry refer to that same field.

In notebooks 2.2 and 2.3, `plot_SourceInfoFromCatalog(..., save_pdf=True)` exports the two-page source report. A PDF path can be passed as `save_pdf="source.pdf"`, and `show=False` suppresses display while saving.

Notebook 2.2 estimates spatial 5-sigma depth from source-masked blank apertures of 3-arcsec diameter on a 10 × 10 grid. Combined depth figures under `output/singlephot/analysis/` select the deepest valid field at each sky position in overlaps. FINALCAT recovery products cover the union of the COSMOS footprints. Blank-aperture sensitivity and the per-band 50%, 80%, and 90% recovery depths are separate measurements.

## Current Status

This project is not complete.

Completed or partially completed:

- Research proposal and scientific motivation.
- Data organization for 7DT COSMOS observations.
- WCS and valid-mask checks.
- SourceXtractor++ configuration for 7DT photometry.
- Per-field COSMOS1–4 single-band photometry and master-catalog workflow.
- FINALCAT union matching and SDSS DR17 Petrosian-photometry integration.
- FINALCAT-centered COSMOS1–4 catalog construction with field measurements, SDSS reconciliation, and explicit FITS mask preservation.
- Band-count/RM diagnostics, per-band recovery curves, and positional-quality products.
- Blank-aperture depth maps and combined-field deepest-coverage maps.
- Shared source inspection with two-page A4 PDF reports.
- Experimental multi-band catalog and diagnostic SED plots.
- Initial CIGALE configuration and 7DT transmission curves.

Still in progress:

- Fully reproducible pipeline scripts.
- Further scientific quality assessment of the combined catalog and selection for SED fitting.
- Final CIGALE SED fitting setup.
- Photo-z estimation.
- Quantitative validation against external redshift catalogs.
- Weak-lensing applicability assessment.

## Git Policy

This repository tracks numbered analysis notebooks, research code, configuration files, filter curves, proposal sources, and input CSVs.

Within `Notebooks/`, only notebooks with a numeric prefix and underscore (for example `1_multibandphot.ipynb` or `2.1_singlephotcatalog_visualization.ipynb`) are versioned. `Proposal.ipynb` and `sketchbook.ipynb` remain available locally but are not tracked. The shared Python helper `Notebooks/sourceinfo_report.py` remains tracked.

The main directories have these roles:

| Directory | Purpose | Git policy |
| --- | --- | --- |
| `Notebooks/` | Numbered analysis workflow and shared source-report helper | Numbered notebooks and code tracked; unnumbered notebooks local |
| `seppconfig/` | SourceXtractor++ photometry configurations | Tracked |
| `sed_fitting/` | CIGALE configuration and filter transmission curves | Tracked |
| `7DTonCOSMOS_tex/Proposal/` | Proposal LaTeX source, bibliography, and source image | Sources tracked; compiled PDF and build files local |
| `Report/` | Report text, input CSVs, and figure/deck builders | Source text/code/inputs tracked; generated figures, PDFs, slides, and build state local |
| `COSMOS_data/` | COSMOS1–4 images and reference catalog | Local only |
| `output/` | Photometric catalogs, matching results, depth grids, and diagnostics | Local only |
| `cigale/` | Local/vendor CIGALE installation | Local only |
| `_tmp/` | Temporary diagnostic products | Local only |

The following are intentionally ignored:

- `COSMOS_data/`: large FITS observations and reference catalogs.
- `output/`: generated catalogs, plots, and intermediate products.
- `cigale/`: local/vendor copy of CIGALE.
- Unnumbered notebooks under `Notebooks/`, including `Proposal.ipynb` and `sketchbook.ipynb`.
- Generated PNG/PDF/PPTX reports, report build/cache directories, `band_hist_depth/` products, and `_tmp/`.
- Python caches, Jupyter checkpoints, OS files, IDE settings, and LaTeX build artifacts.

If large data products need to be versioned later, use a data-management tool such as Git LFS, DVC, or an external archive rather than regular git commits.

## References

This project follows the scientific context of photometric redshift validation and medium-band surveys, including previous 7DT photo-z forecasting work by Ko et al. (2025). Full references are listed in `7DTonCOSMOS_tex/Proposal/references.bib`.
