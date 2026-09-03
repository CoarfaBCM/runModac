# =============================================================================
# splitLipidomics.R
#
# Reusable functions for splitting a targeted-lipidomics export into per-class
# workbooks, grouped by internal standard (ISTD).
#
# This file contains function definitions ONLY - sourcing it has no side effects.
#
#   detectIstdGroups(row1, row2)                  -> data.frame of ISTD groups
#   splitLipidomics(inputFile, outDir, ...)       -> writes workbooks, returns groups
#   writeNormSettings(groups, outDir, ...)        -> writes runModac settings files
#   prepLipidomics(inputFile, splitDir, ...)      -> both writers in one call
#
# The writer functions take a 'methods' argument - "istd", "iqr", or both (the
# default) - so a project can prepare the input workbook and settings file for
# either normalization strategy on its own:
#
#   prepLipidomics(f, splitDir, settingsDir)                      # both
#   prepLipidomics(f, splitDir, settingsDir, methods = "iqr")     # IQR only
#   prepLipidomics(f, splitDir, settingsDir, methods = "istd")    # ISTD only
#
# Input layout assumed (both metadata rows required):
#   row 1  = ISTD roster: for every lipid column, the name of the ISTD used for
#            that lipid's class. The ISTD's own column is typically blank here.
#   row 2  = column names (lipid compounds; the ISTD columns carry the ISTD name)
#   row 3+ = data, one row per sample/injection
#   col 1  = sample name
#
# ISTD detection is derived entirely from the data - no hardcoded column
# positions and no dependence on cell highlighting:
#   1. The roster (unique non-NA values of row 1, after column 1) is the
#      definitive list of ISTDs used in the file.
#   2. Each roster entry is located in row 2 to find that ISTD's own column.
#   3. Sorted ISTD column positions become the group boundaries: an ISTD owns
#      every column after it up to (but not including) the next ISTD.
# =============================================================================


# ---- helpers ----------------------------------------------------------------

# Reduce a name to a comparison key: drop the ISTD marker (in any form and at
# any position), drop all whitespace, lowercase. Everything that identifies the
# compound - parentheses, ':', '_', '-', '+', '.', '/', '|' - is preserved.
istdKey <- function(x) {
  x <- gsub("\\(*ISTD\\)*", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  tolower(x)
}

# Make a string safe for use as an Excel sheet name / R column name.
cleanName <- function(x, maxLen = 31) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  substr(x, 1, maxLen)
}

# Validate the 'methods' argument shared by the writer functions.
# Accepts "istd", "iqr", or both (case-insensitive); errors on anything else.
checkMethods <- function(methods) {
  methods <- unique(tolower(as.character(methods)))
  if (!length(methods)) stop("At least one method must be requested ('istd' and/or 'iqr').")
  bad <- setdiff(methods, c("istd", "iqr"))
  if (length(bad)) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Choose 'istd', 'iqr', or both.")
  }
  methods
}


# ---- ISTD group detection ---------------------------------------------------

# Identify the ISTD groups from the two metadata rows.
#
# row1, row2 : character vectors of equal length (the first two rows of the file)
# returns    : data.frame with one row per ISTD group -
#              istd_col   column index of the ISTD itself
#              istd_name  the ISTD name as written in row 2
#              first,last column range of the lipids belonging to that ISTD
#              n_lipids   number of lipid columns in the group
#
# Errors if a roster entry cannot be matched to exactly one column, or if two
# columns are indistinguishable after normalization. Warns if a match required
# normalization (i.e. row 1 and row 2 spell the ISTD differently) or if a group
# ends up with no lipids.
detectIstdGroups <- function(row1, row2) {

  stopifnot(length(row1) == length(row2))
  n <- length(row1)
  if (n < 3) stop("Need at least 3 columns (sample name + one ISTD + one lipid).")

  keys2 <- istdKey(row2)

  # Guardrail: normalization must not make two different columns identical,
  # otherwise a roster entry could silently resolve to the wrong column.
  nonEmpty <- keys2[!is.na(keys2) & nzchar(keys2)]
  if (anyDuplicated(nonEmpty)) {
    dupKeys <- unique(nonEmpty[duplicated(nonEmpty)])
    stop("Column names are ambiguous after normalization: ",
         paste(dupKeys, collapse = ", "),
         ". Resolve the duplicate names before splitting.")
  }

  # Step 1: the roster - unique non-NA row 1 values after the sample column.
  roster <- unique(row1[2:n])
  roster <- roster[!is.na(roster)]
  if (!length(roster)) {
    stop("No ISTD names found in row 1. Check that row 1 holds the ISTD roster.")
  }

  # Step 2: locate each roster entry in row 2. Exact match is preferred; the
  # normalized fallback covers files where the ISTD marker is written
  # inconsistently between the two rows.
  istdCol <- integer(length(roster))
  for (i in seq_along(roster)) {
    lab <- roster[i]
    hit <- which(row2 == lab)

    if (length(hit) == 0) {
      hit <- which(keys2 == istdKey(lab))
      if (length(hit) == 1) {
        warning("ISTD matched by normalization (row 1 and row 2 differ):\n",
                "  row 1: ", lab, "\n",
                "  row 2: ", row2[hit], call. = FALSE)
      }
    }

    if (length(hit) != 1) {
      stop("ISTD '", lab, "' matched ", length(hit),
           " columns in row 2 (expected exactly 1).",
           if (length(hit) > 1) paste0(" Columns: ", paste(hit, collapse = ", ")) else "")
    }
    istdCol[i] <- hit
  }

  # Step 3: sorted ISTD positions define the group boundaries.
  istdCol <- sort(istdCol)
  groups <- data.frame(
    istd_col  = istdCol,
    istd_name = row2[istdCol],
    first     = istdCol + 1L,
    last      = c(istdCol[-1] - 1L, n),
    stringsAsFactors = FALSE
  )
  groups$n_lipids <- groups$last - groups$first + 1L

  # Validation: every lipid in a group must point back to that group's ISTD.
  for (i in seq_len(nrow(groups))) {
    if (groups$n_lipids[i] < 1) {
      warning("ISTD '", groups$istd_name[i],
              "' (column ", groups$istd_col[i], ") has no lipid columns.",
              call. = FALSE)
      next
    }
    span <- row1[groups$first[i]:groups$last[i]]
    labs <- unique(span[!is.na(span)])
    if (length(labs) != 1) {
      stop("Columns ", groups$first[i], "-", groups$last[i],
           " reference ", length(labs), " different ISTDs in row 1: ",
           paste(labs, collapse = " | "),
           ". Group boundaries are inconsistent.")
    }
    if (istdKey(labs) != istdKey(groups$istd_name[i])) {
      stop("Columns ", groups$first[i], "-", groups$last[i],
           " reference ISTD '", labs, "' but the group's ISTD column is '",
           groups$istd_name[i], "'.")
    }
  }

  groups
}


# ---- splitting --------------------------------------------------------------

# Split a lipidomics export into one sheet per ISTD group.
#
# inputFile      : path to the source .xlsx
# outDir         : directory for the output workbook(s) (created if absent)
# methods        : which workbook(s) to write -
#                    "istd" keeps each group's ISTD column (written last)
#                    "iqr"  drops the ISTD column
#                  default c("istd", "iqr") writes both
# istdFile,
# iqrFile        : output filenames for the two workbooks
# stripDatePrefix: remove a leading date prefix (e.g. "4-13-26-26-") from sample names
# sheet          : sheet index or name to read from inputFile
#
# returns (invisibly) the group table from detectIstdGroups(), with an added
# 'sheet_name' column giving the sheet each group was written to. The table is
# returned regardless of which workbooks were requested, so writeNormSettings()
# can always consume it.
splitLipidomics <- function(inputFile,
                            outDir,
                            methods = c("istd", "iqr"),
                            istdFile = "output_lipidomics_split-ISTD.xlsx",
                            iqrFile = "output_lipidomics_split-IQR.xlsx",
                            stripDatePrefix = TRUE,
                            sheet = 1) {

  methods   <- checkMethods(methods)
  writeIstd <- "istd" %in% methods
  writeIqr  <- "iqr"  %in% methods

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required.")
  }
  if (!file.exists(inputFile)) stop("Input file not found: ", inputFile)
  if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

  raw <- openxlsx::read.xlsx(inputFile, sheet = sheet, colNames = FALSE)
  if (nrow(raw) < 3) stop("Expected at least 3 rows (2 metadata rows + data).")

  row1 <- as.character(unlist(raw[1, ]))
  row2 <- as.character(unlist(raw[2, ]))
  dat  <- raw[3:nrow(raw), , drop = FALSE]

  # read.xlsx() returns every column as character because the two metadata rows
  # are text; restore the measurement columns to numeric so the output workbooks
  # hold real numbers rather than text-formatted ones.
  coerced <- 0
  for (j in 2:ncol(dat)) {
    before <- is.na(dat[[j]])
    dat[[j]] <- suppressWarnings(as.numeric(dat[[j]]))
    coerced <- coerced + sum(is.na(dat[[j]]) & !before)
  }
  if (coerced > 0) {
    warning(coerced, " cell(s) could not be read as numeric and became NA.",
            call. = FALSE)
  }

  groups <- detectIstdGroups(row1, row2)

  message(sprintf("Detected %d ISTD groups covering %d lipid columns (%d samples).",
                  nrow(groups), sum(groups$n_lipids), nrow(dat)))

  # Sheet names come from the ISTD name; de-duplicate if cleaning collides.
  groups$sheet_name <- cleanName(groups$istd_name)
  if (anyDuplicated(groups$sheet_name)) {
    groups$sheet_name <- make.unique(groups$sheet_name, sep = "_")
    groups$sheet_name <- substr(groups$sheet_name, 1, 31)
  }

  # Sample names: strip a leading date prefix, then make them R-safe.
  samples <- as.character(dat[[1]])
  if (stripDatePrefix) samples <- gsub("^\\d+-\\d+-\\d+-\\d+-", "", samples)
  samples <- cleanName(samples, maxLen = 1000L)

  wbIstd <- if (writeIstd) openxlsx::createWorkbook() else NULL
  wbIqr  <- if (writeIqr)  openxlsx::createWorkbook() else NULL

  for (i in seq_len(nrow(groups))) {
    if (groups$n_lipids[i] < 1) next

    lipidCols <- groups$first[i]:groups$last[i]
    sheetName <- groups$sheet_name[i]

    # ISTD column goes last, matching what the normalization settings expect.
    block <- data.frame(dat[, lipidCols, drop = FALSE],
                        dat[, groups$istd_col[i], drop = FALSE],
                        check.names = FALSE, stringsAsFactors = FALSE)
    colnames(block) <- c(row2[lipidCols], groups$istd_name[i])

    if (writeIstd) {
      out <- cbind(Sample.Name = samples, block, stringsAsFactors = FALSE)
      openxlsx::addWorksheet(wbIstd, sheetName)
      openxlsx::writeData(wbIstd, sheetName, out, colNames = TRUE)
    }
    if (writeIqr) {
      out <- cbind(Sample.Name = samples,
                   block[, -ncol(block), drop = FALSE],
                   stringsAsFactors = FALSE)
      openxlsx::addWorksheet(wbIqr, sheetName)
      openxlsx::writeData(wbIqr, sheetName, out, colNames = TRUE)
    }

    message(sprintf("  [%2d] %-31s %3d lipids + 1 ISTD", i, sheetName, groups$n_lipids[i]))
  }

  if (writeIstd) {
    p <- file.path(outDir, istdFile)
    openxlsx::saveWorkbook(wbIstd, p, overwrite = TRUE)
    message("Wrote ", p)
  }
  if (writeIqr) {
    p <- file.path(outDir, iqrFile)
    openxlsx::saveWorkbook(wbIqr, p, overwrite = TRUE)
    message("Wrote ", p)
  }

  invisible(groups)
}


# ---- normalization settings -------------------------------------------------

# Write the runModac "settings" workbooks for ISTD and IQR normalization.
#
# The ISTD settings sheet lists, for each split sheet, the position of that
# sheet's ISTD column. splitLipidomics() puts the ISTD last, so that position is
# n_lipids + 2 (sample-name column + lipids + the ISTD itself).
#
# groups   : the table returned by splitLipidomics()
# outDir   : directory for the settings workbooks (created if absent)
# methods  : which settings file(s) to write - "istd", "iqr", or both (default)
# cvCutoff : value for 'cv_cutoff_internal_standard' in the ISTD settings
writeNormSettings <- function(groups,
                              outDir,
                              methods = c("istd", "iqr"),
                              cvCutoff = 2,
                              istdFile = "settings-ISTD.xlsx",
                              iqrFile = "settings-IQR.xlsx",
                              inputSpace = "linear",
                              diffSpace = "log2") {

  methods <- checkMethods(methods)
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required.")
  }
  if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

  # Row order matters: runModac reads the fixed settings block by position, so
  # the tab -> ISTD_column mapping must begin on row 11 (after 'tab' on row 10).
  base <- function(norm, cv) {
    data.frame(
      key   = c("padj_method", "padj_cutoff", "compute_log2fc", "linear_fc_cutoff",
                "QC_row", "cv_cutoff_internal_standard", "normalization",
                "input_data_space", "differential_analysis_space"),
      value = c(NA, NA, NA, NA, NA, cv, norm, inputSpace, diffSpace),
      stringsAsFactors = FALSE
    )
  }

  writeOne <- function(df, path) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "settings")
    openxlsx::writeData(wb, "settings", df, colNames = FALSE)

    # The value column is character (it mixes text and numbers), so openxlsx
    # would store every cell as text. runModac does arithmetic on the numeric
    # settings directly, so rewrite those cells as real numbers.
    num <- suppressWarnings(as.numeric(df$value))
    for (r in which(!is.na(num))) {
      openxlsx::writeData(wb, "settings", num[r],
                          startRow = r, startCol = 2, colNames = FALSE)
    }

    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    message("Wrote ", path)
  }

  written <- list()

  # IQR: settings block only.
  if ("iqr" %in% methods) {
    iqr <- rbind(base("iqr", NA),
                 data.frame(key = "tab", value = "ISTD_column", stringsAsFactors = FALSE))
    writeOne(iqr, file.path(outDir, iqrFile))
    written$iqr <- file.path(outDir, iqrFile)
  }

  # ISTD: settings block plus one row per sheet giving its ISTD column position.
  if ("istd" %in% methods) {
    if (is.null(groups$sheet_name)) {
      stop("'groups' has no 'sheet_name' column - pass the table returned by splitLipidomics().")
    }
    keep <- groups$n_lipids >= 1
    istd <- rbind(base("istd", cvCutoff),
                  data.frame(key = "tab", value = "ISTD_column", stringsAsFactors = FALSE),
                  data.frame(key   = groups$sheet_name[keep],
                             value = as.character(groups$n_lipids[keep] + 2L),
                             stringsAsFactors = FALSE))
    writeOne(istd, file.path(outDir, istdFile))
    written$istd <- file.path(outDir, istdFile)
  }

  invisible(written)
}


# ---- one-call convenience wrapper -------------------------------------------

# Split the export and write the matching settings file(s) in one call, so the
# requested methods stay in sync between the two steps.
#
# inputFile   : path to the source .xlsx
# splitDir    : directory for the split workbook(s)
# settingsDir : directory for the settings file(s) (defaults to splitDir)
# methods     : "istd", "iqr", or both (default)
#
# returns (invisibly) the group table from splitLipidomics().
prepLipidomics <- function(inputFile,
                           splitDir,
                           settingsDir = splitDir,
                           methods = c("istd", "iqr"),
                           cvCutoff = 2,
                           stripDatePrefix = TRUE,
                           sheet = 1) {

  methods <- checkMethods(methods)
  groups <- splitLipidomics(inputFile = inputFile,
                            outDir = splitDir,
                            methods = methods,
                            stripDatePrefix = stripDatePrefix,
                            sheet = sheet)
  writeNormSettings(groups = groups,
                    outDir = settingsDir,
                    methods = methods,
                    cvCutoff = cvCutoff)
  invisible(groups)
}
