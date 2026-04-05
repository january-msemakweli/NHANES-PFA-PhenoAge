# Renal Clearance and Albumin Binding May Distort Cross-Sectional Associations Between PFAS and Phenotypic Age Acceleration

Reproducible R code and data-extraction pipeline for the analysis of per- and polyfluoroalkyl substance (PFAS) mixtures and Phenotypic Age (PhenoAge) acceleration using NHANES 2005--2018, with a focus on toxicokinetic--algorithm entanglement mediated by chronic kidney disease (CKD) and serum albumin.

## Overview

Cross-sectional environmental epidemiology may be vulnerable to bias when exposure biomarkers are governed by physiologic processes that are also embedded within composite clinical aging algorithms. This project evaluates whether serum PFAS mixture associations with PhenoAge acceleration could reflect toxicokinetic--algorithm overlap related to protein binding and renal clearance.

**Key findings:**

- A one-quantile increase in the four-chemical PFAS mixture (PFOS, PFOA, PFNA, PFHxS) was associated with **lower** PhenoAge acceleration (beta = -1.195; 95% CI: -1.540, -0.850).
- The inverse association was nearly three times stronger among participants with CKD (beta = -2.225) than among those without (beta = -0.773), with a statistically significant interaction (p = 0.015).
- Removing albumin and creatinine from the PhenoAge algorithm attenuated the association (beta = -0.698), consistent with shared physiologic determinants driving the signal rather than intrinsic anti-aging effects.

## Repository Structure

```
.
├── Mixtures_PhenoAge_CKD_NHANES_Analysis.Rmd   # Primary analysis (data extraction, PhenoAge
│                                                 # derivation, qgcomp mixture modeling,
│                                                 # CKD interaction, sensitivity analyses)
├── Participant F-Chart.R                        # Generates the participant flow diagram (Figure 2)
├── full_nhanes_codebook_extraction.R            # Extracts and audits NHANES codebooks across
│                                                 # all cycles for variable harmonization
├── NHANES_CODEBOOK_AUDIT.txt                    # Output of the codebook extraction script
├── DAG 2.png                                    # Directed Acyclic Graph (Figure 1)
└── README.md
```

## Data Source

All data are publicly available from the [National Health and Nutrition Examination Survey (NHANES)](https://www.cdc.gov/nchs/nhanes/index.htm), managed by the CDC's National Center for Health Statistics. Data are downloaded and harmonized at runtime via the [`nhanesA`](https://cran.r-project.org/package=nhanesA) R package. No local data files need to be pre-downloaded.

**Included survey cycles:** 2005--2006, 2007--2008, 2009--2010, 2015--2016, 2017--2018 (cycles lacking C-reactive protein measurements were excluded).

**Final analytic sample:** N = 3,058 adults aged 20 and older with complete PFAS exposure, PhenoAge biomarker, CKD status, survey design, and covariate data.

## Methods

### PhenoAge Derivation

Biological age is estimated using the PhenoAge algorithm (Levine et al., 2018), a Gompertz-based mortality model incorporating chronological age and nine clinical biomarkers: albumin, creatinine, glucose, C-reactive protein, lymphocyte percentage, mean corpuscular volume, red cell distribution width, alkaline phosphatase, and white blood cell count. PhenoAge acceleration is defined as the survey-weighted residual from regressing PhenoAge on chronological age.

### Exposure Assessment

The primary exposure is a four-chemical legacy PFAS mixture (PFOS, PFOA, PFNA, PFHxS). Total PFOS and PFOA are reconstructed from linear and branched isomers for cycles 2015--2018. Concentrations are log2-transformed and ranked into within-cycle quantiles (q = 4) to account for secular declines. The mixture quantile score is the mean of the four within-cycle quantile ranks.

### Statistical Analysis

- **Quantile g-computation** (`qgcomp`) estimates the joint effect of the PFAS mixture on PhenoAge acceleration.
- **Survey-weighted linear regression** (`svyglm`) with CKD interaction terms evaluates effect modification by renal function.
- **Leave-one-out sensitivity analyses** recalculate PhenoAge excluding creatinine and/or albumin to assess toxicokinetic--algorithm overlap.
- All models adjust for sex, race/ethnicity, poverty-income ratio, BMI, smoking status, and NHANES survey cycle.

### CKD Definition

CKD is defined per KDIGO guidelines as eGFR < 60 mL/min/1.73 m² (2021 CKD-EPI equation without race) and/or urine albumin-to-creatinine ratio >= 30 mg/g.

## Requirements

**R version:** >= 4.0

**R packages:**

```r
install.packages(c(
  "nhanesA", "dplyr", "tidyr", "purrr", "stringr",
  "survey", "splines", "ggplot2", "broom", "knitr",
  "qgcomp", "forcats", "tibble"
))
```

## Reproducing the Analysis

1. Clone the repository:

   ```bash
   git clone https://github.com/january-msemakweli/NHANES-PFA-PhenoAge.git
   cd NHANES-PFA-PhenoAge
   ```
2. Install the required R packages (see above).
3. Open and knit `Mixtures_PhenoAge_CKD_NHANES_Analysis.Rmd` in RStudio, or render from the command line:

   ```r
   rmarkdown::render("Mixtures_PhenoAge_CKD_NHANES_Analysis.Rmd")
   ```

   The script downloads all NHANES tables at runtime via the `nhanesA` API. An internet connection is required on first execution. Runtime depends on network speed and NHANES server availability.
4. To regenerate the participant flow chart:

   ```r
   source("Participant F-Chart.R")
   ```
5. To regenerate the codebook audit:

   ```r
   source("full_nhanes_codebook_extraction.R")
   ```

## Survey Design

PFAS is measured in a one-third environmental subsample of NHANES. Cycle-specific subsample weights (e.g., `WTSA2YR`, `WTSB2YR`, `WTSC2YR`) are coalesced into a unified weight and divided by the number of pooled cycles (k = 5) per NCHS analytic guidelines. All variance estimation uses the full complex survey design (strata and PSUs).

## Citation

If you use this code or findings, please cite the accompanying manuscript:

> Msemakweli, JG & Dickerson AS. Renal Clearance and Albumin Binding May Distort Cross-Sectional Associations Between PFAS and Phenotypic Age Acceleration. *[Journal TBD]*.

## Ethics

NHANES has been approved by the NCHS Ethics Review Board. All participants provided written informed consent. IRB documentation is available at [https://www.cdc.gov/nchs/nhanes/about/erb.html](https://www.cdc.gov/nchs/nhanes/about/erb.html).

## License

This project is provided for academic and research purposes. NHANES data are in the public domain.
