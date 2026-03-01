#!/usr/bin/env Rscript

# ==========================================================
# FULL NHANES CODEBOOK EXTRACTION SCRIPT
# Generates: NHANES_CODEBOOK_AUDIT.txt
# ==========================================================

suppressPackageStartupMessages({
  library(nhanesA)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

output_file <- "NHANES_CODEBOOK_AUDIT.txt"
if (file.exists(output_file)) invisible(file.remove(output_file))

cat(
  "NHANES FULL CODEBOOK AUDIT\n",
  "Generated: ", as.character(Sys.time()), "\n",
  "Goal: harmonization-ready metadata by cycle/table/variable.\n\n",
  "NOTE: Analysis starts at 2003-2004. Cycle C uses Lxx_ naming (L40, L25, L06CRP, L06ALB, L10 glucose, L24PFC, L06BMT, L06UHM); later cycles use BIOPRO/CBC/CRP/GLU/PFC/PFAS. CRP in 2015+ is HSCRP.\n\n",
  file = output_file
)

# ----------------------------------------------------------
# 1) Load workspace context (if available)
# ----------------------------------------------------------
if (file.exists(".RData")) {
  load(".RData")
}

safe_get <- function(name, default = NULL) {
  if (exists(name, inherits = TRUE)) get(name, inherits = TRUE) else default
}

default_cycles <- c("2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")  # PFAS mixture: start at 2003-2004 (exclude surplus/pooled 1999-2000, 2001-2002)
config_obj <- safe_get("config", list())
# Use default_cycles so audit always covers all intended cycles (config in .RData may be stale)
cycles <- default_cycles

pfas_vars <- safe_get("pfas_vars", character(0))
blood_metal_vars <- safe_get("blood_metal_vars", character(0))
urine_metal_vars <- safe_get("urine_metal_vars", character(0))
chosen_exposures <- safe_get("chosen_exposures", character(0))
exposure_vars <- safe_get("exposure_vars", character(0))

# ----------------------------------------------------------
# 2) Define all analysis variables used in this project
# ----------------------------------------------------------
analysis_vars <- unique(c(
  # Demographics
  "RIDAGEYR", "RIAGENDR", "RIDRETH1", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
  # Survey design
  "SDMVSTRA", "SDMVPSU", "WTMEC2YR", "WTINT2YR", "WTSAF2YR",
  # Clinical / covariates
  "BMXBMI", "SMQ020", "SMQ040", "ALQ151", "PAQ650", "PAQ665",
  "DIQ010", "BPQ020", "MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F",
  # PhenoAge biomarkers
  "LBXSAL", "LBXSCR", "LBDSCRSI", "LBXGLU", "LBXSGL",
  "LBXCRP", "LBXHSCRP", "LBXLYPCT", "LBXMCV", "LBXMCVSI", "LBXRW", "LBXRDW",
  "LBXAPSI", "LBXSAPSI", "LBXWBCSI", "LBXWBC",
  # CKD markers
  "URXUMA", "URXUCR", "URDUCR", "URDACT",
  # Exposure sets (from workspace + common names)
  pfas_vars, blood_metal_vars, urine_metal_vars, chosen_exposures, exposure_vars,
  "LBXPFOS", "LBXPFOA", "LBXPFNA", "LBXPFHS", "LBXPFDA",
  "LBXBPB", "LBXBCD", "LBXTHG", "URXUAS", "URXUCD", "URXUPB"
))
analysis_vars <- unique(analysis_vars[!is.na(analysis_vars) & nzchar(analysis_vars)])

cycle_suffix <- function(cycle) {
  map <- c(
    "1999-2000" = "A", "2001-2002" = "B", "2003-2004" = "C", "2005-2006" = "D",
    "2007-2008" = "E", "2009-2010" = "F", "2011-2012" = "G", "2013-2014" = "H",
    "2015-2016" = "I", "2017-2018" = "J", "2019-2020" = "K"
  )
  unname(map[cycle])
}

# Cycle-aware candidates: 2003-2004 (C) used Lxx_ naming; later cycles use BIOPRO/CBC/CRP/PFC/PFAS etc.
table_candidates <- list(
  DEMO = c("DEMO"),
  BMX = c("BMX"),
  SMQ = c("SMQ"),
  ALQ = c("ALQ"),
  PAQ = c("PAQ"),
  DIQ = c("DIQ"),
  BPQ = c("BPQ"),
  MCQ = c("MCQ"),
  BIOPRO = c("L40", "BIOPRO"),
  CBC = c("L25", "CBC"),
  CRP = c("L06CRP", "CRP", "HSCRP"),
  ALB_CR = c("L06ALB", "ALB_CR"),
  GLU = c("L10", "GLU"),
  PFAS = c("L24PFC", "PFAS", "PFC", "P_PFAS"),
  BLOOD_METALS = c("L06BMT", "PBCD"),
  URINE_METALS = c("L06UHM", "UHM", "UM")
)

resolve_cycle_table <- function(cycle, candidates) {
  suf <- cycle_suffix(cycle)
  if (is.na(suf) || is.null(suf) || !nzchar(suf)) return(NULL)
  for (base in candidates) {
    to_try <- paste0(base, "_", suf)
    if (!is.character(to_try)) to_try <- as.character(to_try)
    for (tname in to_try) {
      dat <- tryCatch(suppressWarnings(nhanes(tname)), error = function(e) NULL)
      if (!is.null(dat) && nrow(dat) > 0) {
        return(list(table = tname, data = dat))
      }
    }
  }
  NULL
}

cat("Cycles:\n", paste(cycles, collapse = ", "), "\n\n", file = output_file, append = TRUE)
cat("Total tracked analysis variables: ", length(analysis_vars), "\n\n", file = output_file, append = TRUE)

# ----------------------------------------------------------
# 3) Helpers for dictionary parsing and observed summaries
# ----------------------------------------------------------
first_col <- function(df, patterns) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  nms <- names(df)
  idx <- which(vapply(patterns, function(p) any(str_detect(nms, regex(p, ignore_case = TRUE))), logical(1)))
  if (length(idx) == 0) return(NULL)
  pat <- patterns[idx[1]]
  hit <- nms[str_detect(nms, regex(pat, ignore_case = TRUE))]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

extract_var_dict <- function(cb, varname) {
  # nhanesCodebook() usually returns a named list keyed by variable name.
  if (is.list(cb) && !is.data.frame(cb) && !is.null(cb[[varname]])) {
    item <- cb[[varname]]
    label <- NA_character_
    if (is.list(item)) {
      label <- item[["SAS Label:"]]
      if (is.null(label) || length(label) == 0 || is.na(label)) label <- NA_character_
      dict_tbl <- item[[varname]]
      if (inherits(dict_tbl, "data.frame")) {
        code_col <- first_col(dict_tbl, c("^Code or Value$", "^Code.or.Value$"))
        val_desc_col <- first_col(dict_tbl, c("^Value Description$", "^Value.Description$"))
        if (!is.null(code_col) && !is.null(val_desc_col)) {
          value_df <- dict_tbl %>%
            transmute(
              code_or_value = as.character(.data[[code_col]]),
              value_meaning = as.character(.data[[val_desc_col]])
            ) %>%
            filter(!(is.na(code_or_value) | code_or_value == "") | !(is.na(value_meaning) | value_meaning == "")) %>%
            distinct()
          return(list(label = label, value_df = value_df))
        }
      }
    }
    return(list(label = label, value_df = tibble()))
  }

  # Fallback path if cb is a rectangular data frame.
  if (is.null(cb) || !inherits(cb, "data.frame") || nrow(cb) == 0) {
    return(list(label = NA_character_, value_df = tibble()))
  }
  var_col <- first_col(cb, c("^Variable.Name$", "Variable Name", "^Variable$"))
  desc_col <- first_col(cb, c("^Variable.Description$", "Variable Description", "^SAS Label$"))
  code_col <- first_col(cb, c("^Code.or.Value$", "Code or Value", "^Code$"))
  val_desc_col <- first_col(cb, c("^Value.Description$", "Value Description", "^Value$"))
  if (is.null(var_col)) return(list(label = NA_character_, value_df = tibble()))
  cb_var <- cb %>% filter(.data[[var_col]] == varname)
  if (nrow(cb_var) == 0) return(list(label = NA_character_, value_df = tibble()))
  label <- if (!is.null(desc_col)) unique(na.omit(as.character(cb_var[[desc_col]])))[1] else NA_character_
  if (is.na(label) || length(label) == 0) label <- NA_character_
  value_df <- tibble()
  if (!is.null(code_col) && !is.null(val_desc_col)) {
    value_df <- cb_var %>%
      transmute(
        code_or_value = as.character(.data[[code_col]]),
        value_meaning = as.character(.data[[val_desc_col]])
      ) %>%
      filter(!(is.na(code_or_value) | code_or_value == "") | !(is.na(value_meaning) | value_meaning == "")) %>%
      distinct()
  }
  list(label = label, value_df = value_df)
}

format_num <- function(x) {
  if (!is.finite(x)) return("NA")
  format(round(x, 6), trim = TRUE, scientific = FALSE)
}

# ----------------------------------------------------------
# 4) Loop each cycle and each table; audit variables in use
# ----------------------------------------------------------
for (cy in cycles) {
  cat(
    "\n=====================================================\n",
    "CYCLE: ", cy, "\n",
    "=====================================================\n",
    file = output_file, append = TRUE
  )

  cycle_seen_vars <- character(0)
  resolved_tables <- list()
  for (domain in names(table_candidates)) {
    hit <- resolve_cycle_table(cy, table_candidates[[domain]])
    if (!is.null(hit)) resolved_tables[[domain]] <- hit
  }
  if (length(resolved_tables) == 0) {
    cat("No tables resolved for cycle.\n", file = output_file, append = TRUE)
    next
  }

  for (domain in names(resolved_tables)) {
    tb <- resolved_tables[[domain]]$table
    dat <- resolved_tables[[domain]]$data

    vars_in_table <- intersect(names(dat), analysis_vars)
    if (length(vars_in_table) == 0) next

    cycle_seen_vars <- union(cycle_seen_vars, vars_in_table)

    cb <- tryCatch(nhanesCodebook(tb), error = function(e) NULL)

    cat("\nTABLE: ", tb, "\n", file = output_file, append = TRUE)

    for (v in vars_in_table) {
      d <- extract_var_dict(cb, v)

      vec <- dat[[v]]
      vec_chr <- as.character(vec)
      n_total <- length(vec_chr)
      n_missing <- sum(is.na(vec))
      pct_missing <- round(100 * n_missing / max(1, n_total), 2)

      # Observed unique values and frequency
      val_tab <- sort(table(vec_chr, useNA = "ifany"), decreasing = TRUE)
      obs_top <- head(val_tab, 50)
      obs_unique_n <- length(unique(vec_chr))

      # Numeric summary if coercible
      vec_num <- suppressWarnings(as.numeric(vec_chr))
      n_num <- sum(!is.na(vec_num))

      cat(
        "\n---------------------------------------------\n",
        "Variable: ", v, "\n",
        "Label: ", ifelse(is.na(d$label), "<not found in codebook>", d$label), "\n",
        "N total: ", n_total, "\n",
        "N missing: ", n_missing, "\n",
        "Percent missing: ", pct_missing, "%\n",
        "Unique observed values (including NA token): ", obs_unique_n, "\n",
        file = output_file, append = TRUE
      )

      if (n_num > 0) {
        cat(
          "Numeric summary -> Min: ", format_num(min(vec_num, na.rm = TRUE)),
          " | Max: ", format_num(max(vec_num, na.rm = TRUE)),
          " | Mean: ", format_num(mean(vec_num, na.rm = TRUE)),
          " | SD: ", format_num(sd(vec_num, na.rm = TRUE)),
          " | Median: ", format_num(median(vec_num, na.rm = TRUE)), "\n",
          sep = "", file = output_file, append = TRUE
        )
      }

      if (nrow(d$value_df) > 0) {
        cat("\nDictionary codes/meanings:\n", file = output_file, append = TRUE)
        write.table(
          d$value_df, file = output_file, append = TRUE, row.names = FALSE,
          col.names = FALSE, quote = FALSE, sep = "\t"
        )
      } else {
        cat("\nDictionary codes/meanings: <none listed>\n", file = output_file, append = TRUE)
      }

      cat("\nObserved value frequencies (top 50):\n", file = output_file, append = TRUE)
      write.table(
        data.frame(value = names(obs_top), n = as.integer(obs_top)),
        file = output_file, append = TRUE, row.names = FALSE,
        col.names = FALSE, quote = FALSE, sep = "\t"
      )
    }
  }

  missing_in_cycle <- setdiff(analysis_vars, cycle_seen_vars)
  cat(
    "\nCYCLE SUMMARY: analyzed variables found = ", length(cycle_seen_vars),
    " / ", length(analysis_vars), "\n",
    file = output_file, append = TRUE
  )
  if (length(missing_in_cycle) > 0) {
    cat("Variables not found in any downloaded table for this cycle:\n",
        paste(missing_in_cycle, collapse = ", "), "\n", file = output_file, append = TRUE)
  }
}

cat("\n\n================ END OF AUDIT ================\n", file = output_file, append = TRUE)
message("NHANES_CODEBOOK_AUDIT.txt successfully generated.")
