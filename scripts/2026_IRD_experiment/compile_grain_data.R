# Compile instrument summaries; no additional R packages required.
# Run: Rscript scripts/compile_grain_data.R [source_directory] [output_csv]
# Defaults are independent of the current working directory.
args <- commandArgs(trailingOnly = TRUE)
# Ensure accented headers work even when the inherited locale is unavailable.
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_dir <- dirname(dirname(script_path))
source_dir <- if (length(args) >= 1L) args[1] else "D:/dejeanju"
output_file <- if (length(args) >= 2L) args[2] else
  file.path(project_dir, "data", "processed_data", "grain_data.csv")
if (!dir.exists(source_dir)) stop("Source directory does not exist: ", source_dir)

sample_pattern <- "P[0-9]+_? *(mono|mixt|solo|miuxt)[0-9]+$"
folders <- list.dirs(source_dir, recursive = FALSE, full.names = TRUE)
sample_folders <- folders[grepl(sample_pattern, basename(folders))]
ignored <- setdiff(folders, sample_folders)
if (length(ignored)) message("Ignoring non-sample folders: ",
                             paste(basename(ignored), collapse = ", "))
if (!length(sample_folders)) stop("No sample folders found.")

read_sample <- function(folder) {
  sample_id <- regmatches(basename(folder), regexpr(sample_pattern, basename(folder)))
  files <- list.files(folder, pattern = "_METRO_STATISTIQUES\\.dat$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) != 1L) stop(sample_id, ": expected one summary file, found ", length(files))
  # Instrument exports use Windows-1252, tabs, and decimal commas.
  # Only the first two lines describe the whole sample; later tables are grain statistics.
  lines <- readLines(files[1], n = 2L, encoding = "CP1252", warn = FALSE)
  if (length(lines) != 2L) stop(sample_id, ": incomplete summary")
  lines <- iconv(lines, from = "CP1252", to = "UTF-8")
  headers <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  values <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  count_header <- intersect(c("Nb grains", "Nb graines"), headers)
  if (length(count_header) != 1L) stop(sample_id, ": missing or ambiguous grain count")
  wanted <- c("Masse totale mesur\u00e9e", count_header, "PMG")
  positions <- match(wanted, headers)
  if (anyNA(positions) || any(positions > length(values)))
    stop(sample_id, ": missing summary fields")
  numbers <- suppressWarnings(as.numeric(sub(",", ".", values[positions], fixed = TRUE)))
  if (any(!is.finite(numbers)) || any(numbers < 0) || numbers[2] != floor(numbers[2]))
    stop(sample_id, ": invalid mass, grain count, or PMG")
  data.frame(sample_id = sample_id, "Masse totale mesur\u00e9e" = numbers[1],
             "Nb grains" = numbers[2], PMG = numbers[3], check.names = FALSE)
}

result <- do.call(rbind, lapply(sort(sample_folders), read_sample))
if (anyDuplicated(result$sample_id)) {
  # Keep repeats as separate rows, as requested. Preserve folder provenance in outputs.
  review <- result
  review$source_folder <- basename(sort(sample_folders))
  review_file <- file.path(project_dir, "outputs", "grain_measurements_review.csv")
  dir.create(dirname(review_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(review, review_file, row.names = FALSE, fileEncoding = "UTF-8")
  message("Keeping repeated sample IDs as separate rows; provenance saved to ", review_file)
}
result <- result[order(result$sample_id), ]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
# Standard comma-separated CSV, decimal points, UTF-8 column names.
write.csv(result, output_file, row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote ", nrow(result), " measurements (", length(unique(result$sample_id)),
        " distinct sample IDs) to ", output_file)
