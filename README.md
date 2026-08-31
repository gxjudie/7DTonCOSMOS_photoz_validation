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

## Method Summary

The intended analysis workflow is:

1. Inspect calibrated 7DT images and WCS consistency.
2. Build valid masks for regions covered by all usable bands.
3. Generate single-band and/or multi-band detection catalogs using SourceXtractor++.
4. Match catalogs across 7DT filters to construct a master multi-band photometric catalog.
5. Convert measured magnitudes to fluxes and build 20-band SEDs.
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
│   ├── 1_multibandphot.ipynb
│   ├── 2_singphotandmatching.ipynb
│   └── 3_validation.ipynb
├── seppconfig/
│   ├── multiphotconfig.config
│   ├── multiphotconfig.py
│   └── singlephot_cosmos1_m*.config / .py
├── sed_fitting/
│   ├── pcigale.ini
│   ├── pcigale.ini.spec
│   └── transmission_curve/
├── 7DTonCOSMOS_tex/
│   └── Proposal/
├── COSMOS_data/        # local data, ignored by git
├── output/             # generated products, ignored by git
└── cigale/             # local/vendor CIGALE copy, ignored by git
```

## Notebook Roles

- `0_preReaserch.ipynb`: preliminary checks, WCS consistency tests, valid-mask construction, S/N checks, and inspection of the external COSMOS catalog.
- `1_multibandphot.ipynb`: multi-band photometry experiments, chi/chi-mean detection-image construction, SourceXtractor++ execution, and SED plot checks.
- `2_singphotandmatching.ipynb`: single-band photometry and catalog matching workflow.
- `3_validation.ipynb`: reserved for final photo-z validation; currently not completed.

## Current Status

This project is not complete.

Completed or partially completed:

- Research proposal and scientific motivation.
- Data organization for 7DT COSMOS observations.
- WCS and valid-mask checks.
- SourceXtractor++ configuration for 7DT photometry.
- Single-band photometry for COSMOS1 filters.
- Experimental multi-band catalog and diagnostic SED plots.
- Initial CIGALE configuration and 7DT transmission curves.

Still in progress:

- Fully reproducible pipeline scripts.
- Unified multi-field photometric catalog construction.
- Final CIGALE SED fitting setup.
- Photo-z estimation.
- Quantitative validation against external redshift catalogs.
- Weak-lensing applicability assessment.

## Git Policy

This repository tracks research code, notebooks, configuration files, filter curves, and proposal files.

The following are intentionally ignored:

- `COSMOS_data/`: large FITS observations and reference catalogs.
- `output/`: generated catalogs, plots, and intermediate products.
- `cigale/`: local/vendor copy of CIGALE.
- Python caches, Jupyter checkpoints, OS files, IDE settings, and LaTeX build artifacts.

If large data products need to be versioned later, use a data-management tool such as Git LFS, DVC, or an external archive rather than regular git commits.

## References

This project follows the scientific context of photometric redshift validation and medium-band surveys, including previous 7DT photo-z forecasting work by Ko et al. (2025). Full references are listed in `7DTonCOSMOS_tex/Proposal/references.bib`.
