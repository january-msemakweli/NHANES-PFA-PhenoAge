# Handling of Variables: Source, Recoding, and Methodological Notes

**Project:** Mixtures (PFAS) - PhenoAge - CKD analysis using NHANES 2003-2018**Files this document covers:**

- `Mixtures_PhenoAge_CKD_NHANES_Analysis.Rmd` (v1, primary analysis)
- `Mixtures_PhenoAge_CKD_NHANES_Analysis_TableOne_v2.Rmd` (v2, two-part Table 1)
- `Mixtures_PhenoAge_CKD_NHANES_Analysis_TableOne_v3_export.Rmd` (v3, two-part Table 1 + analytic-sample CSV exports)

All three files share the same upstream variable construction (helpers + `{r construct-vars}` chunk + biomarker block). They diverge only in (a) Table 1 layout (v1 vs v2/v3) and (b) the final CSV export step (v3 only). Unless otherwise noted, code references in this document refer to the **v3 export file** (which is the superset).

---

## 1. Conventions

- **Cycles included:** NHANES C-J = 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016, 2017-2018 (8 cycles, 16 calendar years).
  - 1999-2000 (A) and 2001-2002 (B) are intentionally excluded: they lack the standard PFAS panel and use surplus/prototype assays (`SSPFC_A`) that are not weight-compatible with later cycles.
  - See `config$cycles` in `{r config}`.
- **Adult restriction:** `age >= 20` (NHANES adult definition for clinical labs and PFAS panels). Set in `config$min_age` and applied at the start of `{r analytic-sample}`.
- **Source-table joins:** For each cycle, every NHANES table needed for the study is downloaded via `nhanesA::nhanes()` and joined on `SEQN` (`download_cycle()` in `{r data-download}`). All cycles are then column-harmonized (missing columns filled with `NA`) and stacked with `bind_rows` into `nhanes_raw`.
- **Cycle-dependent column guard:** A defensive block at the top of `{r construct-vars}` ensures every cycle-dependent column (e.g. `ALQ142`, `PAD200`) exists in `nhanes_raw` as `NA` if no included cycle had it. This makes the harmonized derivations below robust to changes in `config$cycles`.
- **Coalescing across renamed NHANES variables:** Many biomarkers were renamed across cycles (e.g. `LBXRW` -> `LBXRDW`). The helper `coalesce_num(df, c(...))` returns the first non-`NA` numeric value across the candidate column names per row.
- **Yes/No parsing:** Helpers in `{r helpers}` (`parse_yes_no`, `yesno01`, `parse_smq040`, `parse_smq020`, `parse_alq_ever`, `parse_alq_current_yr`) all map `1 -> 1` (Yes), `2 -> 0` (No), and `7/9` ("Refused"/"Don't know") -> `NA`. They also accept legacy character labels ("Yes"/"No") for cycles where `nhanesA` returns labels.
- **Standardized missing-data semantics:** `NA` always means "information is missing" (not "No"). Refused/Don't know is treated as missing. Skip-pattern responses are recovered to the appropriate non-missing value whenever the gate question can be used to disambiguate (see Smoking, Alcohol).

---

## 2. Identifiers and Cycle Metadata

| Variable  | NHANES source                   | Type | Notes                                                                                         |
| --------- | ------------------------------- | ---- | --------------------------------------------------------------------------------------------- |
| `SEQN`  | `SEQN`                        | int  | NHANES sequence number; unique participant ID within a cycle.                                 |
| `cycle` | derived in `download_cycle()` | char | "2003-2004" through "2017-2018"; used for weight pooling, audits, and the cycle fixed effect. |

---

## 3. Exposures: PFAS panel

**Primary 4-PFAS set** (cycle-stable across 2003-2018, included in every retained cycle):

| Derived variable | NHANES source                                                             | Description                             |
| ---------------- | ------------------------------------------------------------------------- | --------------------------------------- |
| `LBXPFOS`      | `LBXPFOS` (or `LBXMPAH`-suffixed in J cycle, see harmonization notes) | Perfluorooctane sulfonate (ng/mL serum) |
| `LBXPFOA`      | `LBXPFOA`                                                               | Perfluorooctanoate                      |
| `LBXPFNA`      | `LBXPFNA`                                                               | Perfluorononanoate                      |
| `LBXPFHS`      | `LBXPFHS`                                                               | Perfluorohexane sulfonate               |

**Cycle-specific harmonization (J cycle, 2017-2018):** PFOS and PFOA were measured as linear + branched isomer pairs in J. The download routine (`{r data-download}`) sums the linear + branched columns and writes the result back into the legacy `LBXPFOS` / `LBXPFOA` columns so downstream code is cycle-agnostic. See lines around 880-900 of v3.

**Limit-of-detection (LOD) handling:** `apply_lod_rule(df, varname)` in `{r helpers}` (lines ~311-340) implements the standard NHANES PFAS rule:

- If a `<varname>LC` (or `LBD<base>LC`/`URD<base>LC`) detection-limit flag exists and equals 1 (below LOD), the concentration is replaced with `LOD / sqrt(2)`.
- For NHANES tables that store the LOD value directly in the concentration column when LC=1, the rule reduces to dividing the stored value by `sqrt(2)`.
- If the value is still `NA` or non-positive after that, the smallest observed positive LOD/sqrt(2) for that analyte is used as a fallback.
- A per-cycle LOD diagnostic table is printed at the end of `{r exposure-construct}` so reviewers can verify the imputation behavior.

**Log transform & mixture quantization:**

- After LOD substitution, every exposure is transformed as `<exposure>_log2 = log2(value)`.
- For the qgcomp mixture model, each `_log2` exposure is quantized into `q = 4` quartiles (`config$q`). The `mixture_qscore` is the row-mean of these quartile ranks across the 4 PFAS, with rows requiring at least `min_nonmissing_exposures_floor = 2` non-missing exposures.
- See `{r qgcomp-prep}` block (lines ~1220-1360 in v3).

---

## 4. Outcome: PhenoAge and PhenoAge Acceleration

**Construct:** Levine et al. (2018) Phenotypic Age, computed from 9 clinical biomarkers + chronological age via `compute_phenoage_components()` in `{r helpers}` (lines ~241-271).

**Inputs (all coalesced and unit-converted in `{r construct-vars}` lines ~600-609):**

| PhenoAge component      | Variable | NHANES source columns (coalesced)                 | Levine unit   | NHANES unit (after conversion)    |
| ----------------------- | -------- | ------------------------------------------------- | ------------- | --------------------------------- |
| Albumin                 | `alb`  | `LBXSAL`                                        | g/L           | g/dL -> ×10                      |
| Creatinine              | `scr`  | `LBXSCR`, `LBDSCRSI`                          | µmol/L       | mg/dL -> ×88.4                   |
| Glucose (fasting)       | `glu`  | `LBXGLU`, `LBXSGL` (primary uses MEC fasting) | mmol/L        | mg/dL -> ÷18.018                 |
| C-reactive protein      | `crp`  | `LBXCRP`, `LBXHSCRP`                          | mg/L          | mg/L (log of `pmax(crp, 1e-6)`) |
| Lymphocyte %            | `lym`  | `LBXLYPCT`                                      | %             | %                                 |
| Mean corpuscular volume | `mcv`  | `LBXMCV`, `LBXMCVSI`                          | fL            | fL                                |
| RDW                     | `rdw`  | `LBXRW`, `LBXRDW`                             | %             | %                                 |
| Alkaline phosphatase    | `alp`  | `LBXAPSI`, `LBXSAPSI`                         | U/L           | U/L                               |
| WBC                     | `wbc`  | `LBXWBCSI`, `LBXWBC`                          | 10^3 cells/uL | 10^3 cells/uL                     |
| Age                     | `age`  | `RIDAGEYR`                                      | years         | years                             |

**Numeric stability:**

- The linear predictor is clipped to `[-50, 50]`.
- `mort` is clipped to `[1e-12, 1 - 1e-12]` before back-transformation.

**Sensitivity variants:**

- `phenoage_nokid` = PhenoAge with the creatinine term zeroed out (`include_scr = FALSE`). Used to detect kidney-driven biomarker confounding.
- `phenoage_nokid_both` = PhenoAge with both creatinine and albumin terms zeroed out.
- See `{r sensitivity}` block.

**PhenoAge Acceleration (`phenoage_accel`):**

- Residual of `phenoage ~ age` from a survey-weighted (PSU/strata/`weight_combined`) Gaussian linear model fit on the analytic sample.
- A second, unweighted residualization (`phenoage_accel_unweighted`) is computed for sensitivity comparisons.
- See `{r residualize-phenoage}` block (lines ~1056-1085).

---

## 5. CKD: eGFR, ACR, and the binary outcome

**eGFR (`egfr`):** CKD-EPI 2021 race-free equation. Implemented in `compute_egfr_2021()` in `{r helpers}` (lines ~220-235). Inputs: serum creatinine (`scr`, mg/dL), age (years), sex (`RIAGENDR` raw 1=Male / 2=Female).

**ACR (`acr`):** `compute_acr(ualb, ucr, acr_direct)` in `{r helpers}` (lines ~236-240). Logic:

1. If `URDACT` (NHANES-derived ACR, mg/g) is present, use it directly.
2. Otherwise compute as `(ualb [mg/L] / ucr [mg/dL]) * 100` from `URXUMA` and `URXUCR`/`URDUCR`.

**CKD components (in `{r construct-vars}` lines ~616-624):**

| Variable      | Definition                                                                                                                                                                            |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ckd_egfr`  | `1` if `egfr < 60`, else `0`; `NA` if eGFR missing.                                                                                                                           |
| `ckd_acr`   | `1` if `acr >= 30 mg/g`, else `0`; `NA` if ACR missing.                                                                                                                       |
| `ckd`       | `1` if either marker is `1`, `0` if either is observed and neither is `1`, else `NA`. Selectable via `config$ckd_def` ("egfr_or_acr" [primary], "egfr_only", "acr_only"). |
| `ckd_stage` | KDIGO G-stages from `cut(egfr, breaks = c(-Inf,15,30,45,60,90,Inf), labels = c("G5","G4","G3b","G3a","G2","G1"))`, left-closed intervals.                                           |

**Sensitivity (`ckd_sens`):** Same definition with `config$ckd_def_sens` (typically "acr_only" or stricter eGFR threshold). See `{r sensitivity}`.

---

## 6. Demographics

### 6.1 Age (`age`)

- **Source:** `RIDAGEYR` (years).
- **Recoding:** `as.numeric(as.character(RIDAGEYR))` (NHANES top-codes 80+ as 80 in some cycles).
- **Inclusion:** `age >= 20`.

### 6.2 Sex (`sex`)

- **Source:** `RIAGENDR` (1=Male, 2=Female).
- **Recoding:** Mapped to character "Male"/"Female". String fallback parses raw "MALE"/"FEMALE" labels for cycles where `nhanesA` returns labels.

### 6.3 Race/Ethnicity (`race_eth`)

- **Sources:** `RIDRETH3` (preferred, available 2011+ and 2003-2006 with Asian category) coalesced to `RIDRETH1` (legacy 4-level) for older cycles.
- **5-level harmonization** (chosen to match published NHANES descriptive practice and to expose Non-Hispanic Asian as its own row in Table 1):
  - `RIDRETH3 ∈ {1,2}` = "Hispanic" (Mexican-American + Other Hispanic)
  - `RIDRETH3 == 3` = "Non-Hispanic White"
  - `RIDRETH3 == 4` = "Non-Hispanic Black"
  - `RIDRETH3 == 6` = "Non-Hispanic Asian"
  - `RIDRETH3 ∈ {5,7}` = "Other/Multi"
- **Fallback for early cycles** (no RIDRETH3): `RIDRETH1` is used and coded into the same 5 categories, with NHA folded into "Other/Multi" for cycles before the Asian oversample (2003-2010).
- A per-cycle audit table (`Race audit: availability of RIDRETH3 by cycle`) is printed in `{r construct-vars}` to flag cycles relying on the fallback.

### 6.4 Education (`education`)

- **Source:** `DMDEDUC2` (asked of adults 20+).
- **3-level recode:**
  - `1, 2` -> "<HS" (Less than 9th + 9-11th grade)
  - `3` -> "HS/GED"
  - `4, 5` -> ">HS" (Some college + College graduate)
  - `7, 9` (Refused/DK) -> `NA`

### 6.5 Poverty-Income Ratio (`pir`)

- **Source:** `INDFMPIR` (continuous family poverty-income ratio).
- **Recoding:** Pure numeric coercion. NHANES top-codes at 5.0.

---

## 7. Anthropometry

### 7.1 BMI (`bmi`)

- **Source:** `BMXBMI` (kg/m^2, measured in MEC).
- **Recoding:** `as.numeric(as.character(BMXBMI))`.

---

## 8. Behavioral Covariates

### 8.1 Smoking (`smoker3`, `smoker`)

**The problem:** NHANES smoking is a 2-stage instrument:

- `SMQ020` ("Smoked >=100 cigarettes in life?") is the GATE question.
- `SMQ040` ("Do you now smoke cigarettes?") is asked **only** of `SMQ020 == Yes` respondents.

A naive derivation that uses only `SMQ040` (the previous version of this code) silently treats every never-smoker (`SMQ020 == No`, who is never asked `SMQ040`) as `NA` and drops them from any complete-case modeling sample. This is a substantial bias.

**Fix (current implementation, `{r construct-vars}` lines ~512-535):**

```r
smoker3 = {
  ever <- parse_smq020(SMQ020)          # 1=ever-smoker, 0=never-smoker
  curr <- <recode SMQ040 to 0/1>        # 1=current, 0=former; NA if SMQ020==No
  case_when(
    ever == 0 ~ "Never",
    ever == 1 & curr == 1 ~ "Current",
    ever == 1 & curr == 0 ~ "Former",
    TRUE ~ NA_character_
  )
},
smoker = case_when(
  smoker3 == "Current" ~ 1,
  smoker3 %in% c("Never", "Former") ~ 0,
  TRUE ~ NA_real_
)
```

| Derived     | Levels                                        | Used by                                                      |
| ----------- | --------------------------------------------- | ------------------------------------------------------------ |
| `smoker3` | "Never" / "Former" / "Current"                | Table 1 (descriptive).                                       |
| `smoker`  | 0 = Not current (Never + Former), 1 = Current | All adjusted models (qgcomp, sensitivity, marginal effects). |

**Why a binary `smoker` in the regression model:**

1. Comparability to the published literature on PFAS-PhenoAge (Wang et al., Cohen et al., Calafat et al.) which uses Current vs Not-current.
2. Identifiability in a CKD-stratified subgroup (n ≈ 585) where 3-level smoking yields fragile coefficients.
3. The 3-level distribution is preserved descriptively in `smoker3` for Table 1 transparency.
4. Sensitivity analysis with `smoker3` as a 3-level fixed effect is recommended (and easy because both vars are exported).

### 8.2 Alcohol (`alcohol3`, `alcohol`)

**The problem:** The previous derivation used only `ALQ151` ("Ever 4/5+ drinks every day in any 12-month period?"), which:

- Measures *heavy episodic drinking history*, not current drinker status (mismatch with PI's intended construct).
- Has a hard skip pattern: never-drinkers (gated out by `ALQ101`/`ALQ110`) are recorded as `NA`, not 0.
- Was discontinued in the 2017-2018 (J) cycle, dropping the entire J cycle to `NA`.

**Fix (current implementation, `{r construct-vars}`):**

The derivation now uses two cycle-harmonized helpers (`{r helpers}`):

```r
parse_alq_ever(ALQ101, ALQ110, ALQ111)        # lifetime ever-drinker
parse_alq_current_yr(ALQ121, ALQ120Q)         # drank in past 12 months
```

| NHANES question                                   | Variable    | Cycles    | Values                                             |
| ------------------------------------------------- | ----------- | --------- | -------------------------------------------------- |
| ">=12 drinks in any one year ever?"               | `ALQ101`  | 1999-2016 | 1=Yes / 2=No / 7,9=Ref/DK                          |
| ">=12 drinks in lifetime?" (asked when ALQ101=No) | `ALQ110`  | 1999-2016 | 1=Yes / 2=No / 7,9=Ref/DK                          |
| "Ever had a drink of any kind of alcohol?"        | `ALQ111`  | 2017+     | 1=Yes / 2=No / 7,9=Ref/DK                          |
| "Past 12 mo, how often drink alcohol?"            | `ALQ121`  | 2013-2018 | 0=Never in past yr; 1-10=frequencies; 77,99=Ref/DK |
| "How often drink alcohol?" + unit `ALQ120U`     | `ALQ120Q` | 1999-2016 | 0-365 days; 777,999=Ref/DK                         |

`alcohol3` is then defined identically to `smoker3`:

- `Never` if ever-drinker is `0`
- `Current` if ever-drinker `1` AND past-12-mo drinker `1`
- `Former` if ever-drinker `1` AND past-12-mo drinker `0`
- `NA` otherwise (Refused/DK on the gate)

`alcohol` (binary, model covariate) = `1` if `alcohol3 == "Current"`, `0` otherwise.

| Derived      | Levels                                                        | Used by                |
| ------------ | ------------------------------------------------------------- | ---------------------- |
| `alcohol3` | "Never" / "Former" / "Current"                                | Table 1 (descriptive). |
| `alcohol`  | 0 = Not current drinker (Never + Former), 1 = Current drinker | All adjusted models.   |

**Effect on missingness:** Drops from ~70% missing (ALQ151-only definition) to <5% missing (refused/DK only).

### 8.3 Physical Activity (`phys_active`)

**The problem:** NHANES changed the physical activity questionnaire in 2007:

- 1999-2006 (cycles A-D) used `PAD200` (vigorous past 30 days), `PAD320` (moderate past 30 days).
- 2007+ (cycles E onward) replaced these with the GPAQ items `PAQ650` (vigorous recreational), `PAQ665` (moderate recreational).

A `PAQ650/PAQ665`-only derivation drops cycles C and D (2003-2006) entirely.

**Fix (current implementation, `{r construct-vars}`):**

```r
phys_active = case_when(
  parse_yes_no(PAQ650) == 1 | parse_yes_no(PAQ665) == 1 ~ 1,   # GPAQ Yes
  parse_yes_no(PAQ650) == 0 & parse_yes_no(PAQ665) == 0 ~ 0,   # GPAQ No
  parse_yes_no(PAD200) == 1 | parse_yes_no(PAD320) == 1 ~ 1,   # legacy Yes
  parse_yes_no(PAD200) == 0 & parse_yes_no(PAD320) == 0 ~ 0,   # legacy No
  TRUE ~ NA_real_
)
```

The harmonized binary captures *any moderate-or-vigorous physical activity* across either questionnaire format, which is the standard cross-cycle approach in the NHANES literature (e.g. Wang et al. JAMA 2018; Wallace et al. Environ Int 2023).

**Caveat:** The two questionnaires are not strictly identical (GPAQ asks about *recreational* activity; the legacy module asked about activity in the past 30 days regardless of context). The harmonized binary trades some construct precision for cross-cycle comparability. A cycle fixed effect in every model absorbs residual cycle-instrument heterogeneity.

---

## 9. Clinical Covariates

### 9.1 Diabetes (`diabetes`)

- **Source:** `DIQ010` ("Doctor told you have diabetes?"). Codes: `1`=Yes, `2`=No, `3`=Borderline, `7,9`=Refused/DK.
- **Recoding:** `yesno01(DIQ010)` -> `1`/`0`/`NA`. Borderline (`3`) maps to `NA`.
- **Caveat documented in code:** Borderline-diabetics are dropped from the analytic sample. Alternative codings (Borderline -> 1 [any diabetes] or Borderline -> 0 [definite-only]) can be selected if PI prefers; this would shift the analytic `n` up by ~150-200 across cycles.

### 9.2 Hypertension (`hypertension`)

- **Source:** `BPQ020` ("Doctor told you have HBP?"). Codes: `1`=Yes, `2`=No, `7,9`=Refused/DK.
- **Recoding:** `yesno01(BPQ020)` -> `1`/`0`/`NA`. No gate question.
- **Note:** This is *self-reported diagnosis*, not measured BP from MEC. A measured-BP sensitivity definition (`BPXSY1`/`BPXDI1` >= 130/80) can be added if needed.

### 9.3 Cardiovascular Disease history (`cvd`)

- **Source:** Composite of four MCQ items, each "Doctor ever told you had ...":
  - `MCQ160B` = Congestive heart failure
  - `MCQ160C` = Coronary heart disease
  - `MCQ160E` = Heart attack (myocardial infarction)
  - `MCQ160F` = Stroke
- **Recoding:** `1` if ANY component is Yes; `0` if all four are observed and No; `NA` if all four are missing/refused.
- **Caveat:** No gate question among the four. The compound `NA`-handling means a participant with three "No"s and one "Don't know" is still coded `0`, which is the standard NHANES CVD-history convention.

---

## 10. Survey Design

### 10.1 PSU and Strata

- `psu` <- `SDMVPSU`
- `strata` <- `SDMVSTRA`
- Set in `{r construct-vars}` (lines ~630-631) via `coalesce_num`.
- `options(survey.lonely.psu = "adjust")` is set globally.

### 10.2 Weight selection (`weight_2yr`, `weight_combined`)

PFAS were measured in **rotating subsamples** with subsample-specific 2-year weights that change name across cycles:

| Cycle           | PFAS subsample weight                          |
| --------------- | ---------------------------------------------- |
| 2003-2008 (C-E) | `WTSC2YR` (initial PFAS subsample)           |
| 2009-2014 (F-H) | `WTSB2YR` (PFAS subsample B)                 |
| 2015-2018 (I-J) | `WTSA2YR` (PFAS subsample A) or `WTSAF2YR` |

**Harmonization (`{r analytic-sample}` lines ~977-981):** A unified `WTPFAS2YR` is built by coalescing across the four candidate names. The auto-weight resolver (lines ~1003-1015) audits which `WT*` variable has the highest non-missing coverage among participants with the complete primary-PFAS exposure set, then automatically selects it. Generic `WTINT2YR` and `WTMEC2YR` are explicitly excluded from PFAS auto-selection because they over-cover the analytic sample.

**Multi-cycle pooling (`repool_weights()`, line ~1033):** Per NCHS multi-cycle weighting guidance (Mirel et al., 2013; NHANES Analytic Guidelines), pooled-cycle weights must be divided by the number of contributing cycles. `weight_combined = weight_2yr / ncy`, where `ncy` is the number of distinct cycles in the *current* sample (so the divisor changes between `analytic`, `analytic_qg`, `qdat`, etc., as exclusions tighten).

**Survey design objects:**

- `svy_resid` (PhenoAge residualization): `svydesign(id=~psu, strata=~strata, weights=~weight_combined, nest=TRUE, data=analytic)`
- `qdat_svy` (qgcomp interaction model): same template on the qgcomp-eligible subset.
- See `{r residualize-phenoage}`, `{r mixture-models}`, and `{r ckd-stratified}` chunks.

### 10.3 Lonely-PSU policy

Set globally via `options(survey.lonely.psu = "adjust")` in `{r packages}`. Recovery design objects fall back to `svydesign(id=~1, strata=NULL, weights=~weight_combined, ...)` (effectively design-unaware) only when the stratified design errors out on a too-small subgroup; this fallback is logged.

---

## 11. Analytic Sample Construction

The pipeline produces several nested samples:

| Object          | Inclusion criteria                                                                                                                                                  | Approx. n (8 cycles)                       |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------ |
| `nhanes_raw`  | All adults+children downloaded for any included cycle                                                                                                               | ~80,000                                    |
| `dat`         | Same, after variable derivation (no row filtering)                                                                                                                  | ~80,000                                    |
| `analytic`    | +`age >= 20`, non-`NA` `phenoage`, `ckd`, `weight_2yr`, `strata`, `psu`                                                                               | ~15,200                                    |
| `analytic_qg` | + non-`NA` for the 4 primary PFAS (after LOD imputation) and the 6 model covariates (`age`, `sex`, `race_eth`, `education`, `pir`, `bmi`, `smoker`) | ~6,800 (was ~3,058 before the smoking fix) |
| `qdat`        | `analytic_qg` + valid `mixture_qscore`                                                                                                                          | same as above                              |

**Inclusion flag (v3 export):** The full-sample CSV (`analytic_dataset_full_with_inclusion.csv`) carries a column `inclusion_status` = `"Included"` (in `analytic_qg`) or `"Excluded"` (in `analytic` only). This makes Table 1 Part A vs Part B reproducible from the CSV alone.

---

## 12. Variable Revision History (relevant to this project)

| Date            | Variable(s)                                  | Change                                                                                                                                      | Reason                                                                                                                                                                       |
| --------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Initial         | `smoker`                                   | Coded from `SMQ040` only as Current vs Non-current.                                                                                       | Standard short-form NHANES recode; missed gate-question logic.                                                                                                               |
| Smoking fix     | `smoker3` (new), `smoker` (re-derived)   | 3-level Never/Former/Current via `SMQ020 + SMQ040`; binary now Current vs Not-current (Never+Former).                                     | Recovers never-smokers (~3,800 participants who were silently dropped from the complete-case sample). PI-requested justification of binary form is in this doc, section 8.1. |
| p-value display | Figure 2 (forest plot)                       | `geom_text` label now renders `p<0.001` instead of `p=0.000` for very small p-values.                                                 | Display correctness; small p-values were misleadingly shown as zero.                                                                                                         |
| Alcohol fix     | `alcohol3` (new), `alcohol` (re-derived) | 3-level Never/Former/Current via cycle-harmonized `ALQ101/ALQ110/ALQ111` + `ALQ121/ALQ120Q`; binary now Current vs Not-current drinker. | Eliminates the `ALQ151` skip-pattern bias and the J-cycle missingness; aligns with NHANES literature convention.                                                           |
| PA fix          | `phys_active`                              | Added `PAD200/PAD320` legacy fallback for 1999-2006 cycles.                                                                               | Recovers the 2003-2006 cycles that the GPAQ-only definition was dropping.                                                                                                    |

---

## 13. CSV Exports (v3 only)

`{r export-data}` writes two files at the end of v3:

1. `analytic_dataset_complete_case_qgcomp.csv` -- the `analytic_qg` complete-case sample used by qgcomp models. One row per participant, all categorical variables labeled.
2. `analytic_dataset_full_with_inclusion.csv` -- the union of `analytic` (all adults with PhenoAge + CKD + weights) plus an `inclusion_status` flag distinguishing the qgcomp complete-case from the rest.

Both exports use `na = ""` (empty string for missing) and `row.names = FALSE`. Categorical labels are applied by `label_dataset()` (lines ~2380-2400 of v3) so consumers receive human-readable factors rather than 0/1 codes.

---

## 14. Quick-Reference Variable Index

| Domain              | Final variable(s)                                                                          | Source columns                                                                      | Code chunk                                      |
| ------------------- | ------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- | ----------------------------------------------- |
| ID                  | `SEQN`, `cycle`                                                                        | `SEQN`; cycle string                                                              | `{r data-download}`                           |
| Demographics        | `age`, `sex`, `race_eth`, `education`, `pir`                                     | `RIDAGEYR`, `RIAGENDR`, `RIDRETH3`/`RIDRETH1`, `DMDEDUC2`, `INDFMPIR`   | `{r construct-vars}`                          |
| Anthropometry       | `bmi`                                                                                    | `BMXBMI`                                                                          | `{r construct-vars}`                          |
| Smoking             | `smoker3`, `smoker`                                                                    | `SMQ020`, `SMQ040`                                                              | `{r construct-vars}`                          |
| Alcohol             | `alcohol3`, `alcohol`                                                                  | `ALQ101`, `ALQ110`, `ALQ111`, `ALQ120Q`, `ALQ121`                         | `{r construct-vars}`                          |
| Physical activity   | `phys_active`                                                                            | `PAQ650`, `PAQ665`, `PAD200`, `PAD320`                                      | `{r construct-vars}`                          |
| Clinical            | `diabetes`, `hypertension`, `cvd`                                                    | `DIQ010`, `BPQ020`, `MCQ160B/C/E/F`                                           | `{r construct-vars}`                          |
| PhenoAge biomarkers | `alb`, `scr`, `glu`, `crp`, `lym`, `mcv`, `rdw`, `alp`, `wbc`            | See §4 table                                                                       | `{r construct-vars}` lines ~600-609           |
| Outcome             | `phenoage`, `phenoage_accel`, `phenoage_accel_unweighted`                            | derived from biomarkers +`age`                                                    | `{r residualize-phenoage}`                    |
| CKD                 | `ckd_egfr`, `ckd_acr`, `ckd`, `ckd_stage`, `ckd_sens`                            | `LBXSCR`/`LBDSCRSI`, `URXUMA`, `URXUCR`/`URDUCR`, `URDACT`              | `{r construct-vars}` lines ~611-628           |
| Survey design       | `psu`, `strata`, `weight_2yr`, `weight_combined`                                   | `SDMVPSU`, `SDMVSTRA`, `WTSC2YR/WTSB2YR/WTSA2YR/WTSAF2YR`                     | `{r analytic-sample}`                         |
| Exposures           | `LBXPFOS_log2`, `LBXPFOA_log2`, `LBXPFNA_log2`, `LBXPFHS_log2`, `mixture_qscore` | `LBXPFOS`, `LBXPFOA`, `LBXPFNA`, `LBXPFHS` (+ J-cycle isomer harmonization) | `{r exposure-construct}`, `{r qgcomp-prep}` |

---

## 15. Methodological References

- **PhenoAge formula:** Levine ME, Lu AT, Quach A, et al. *An epigenetic biomarker of aging for lifespan and healthspan.* Aging (Albany NY). 2018;10(4):573-591.
- **eGFR equation:** Inker LA, Eneanya ND, Coresh J, et al. *New Creatinine- and Cystatin C-Based Equations to Estimate GFR without Race.* N Engl J Med. 2021;385:1737-1749.
- **NHANES multi-cycle pooling:** Mirel LB, Mohadjer LK, Dohrmann SM, et al. *National Health and Nutrition Examination Survey: Estimation procedures, 2007-2010.* Vital Health Stat 2. 2013;(159):1-17.
- **PFAS LOD substitution (LOD/sqrt(2)):** Hornung RW, Reed LD. *Estimation of average concentration in the presence of nondetectable values.* Appl Occup Environ Hyg. 1990;5(1):46-51.
- **Quantile g-computation (qgcomp):** Keil AP, Buckley JP, O'Brien KM, et al. *A quantile-based g-computation approach to addressing the effects of exposure mixtures.* Environ Health Perspect. 2020;128(4):047004.
