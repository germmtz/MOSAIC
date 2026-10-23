# Run: Rscript scripts/build_trait_files.R
# Requires readxl. Repeated measurements remain separate; source order is preserved.
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl is required.")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
raw <- file.path(root, "data", "raw_data")
processed <- file.path(root, "data", "processed_data")
out <- file.path(root, "outputs")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# French Excel: semicolon fields, decimal commas, UTF-8 BOM for accented comments.
write_output <- function(x, path) {
  for (n in names(x)) if (is.character(x[[n]]))
    x[[n]][!is.na(x[[n]]) & !nzchar(trimws(x[[n]]))] <- NA_character_
  lines <- character()
  text_con <- textConnection("lines", "w", local = TRUE)
  write.csv2(x, text_con, row.names = FALSE, na = "NA")
  close(text_con)
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\r\n"), "\r\n"))), con)
}

reference <- as.data.frame(readxl::read_excel(file.path(raw, "identifiants_echantillons_MOSAIC.xlsx")))
# In the source legend, r = resown; 1 = germinated since 06/02, not resown.
# A blank is treated as no recorded reseeding. Flags apply to pot/genotype groups:
# repeated measurements cannot be assigned to individual plants in mono pots.
reseeding_codes <- tolower(trimws(as.character(reference[[5]])))
if (any(!is.na(reseeding_codes) & !reseeding_codes %in% c("", "r", "1", "0")))
  stop("Unrecognized reseeding code in the reference workbook")
# Carry pot metadata down within each pot, but never carry genotype or notes down.
for (i in seq_len(nrow(reference))[-1]) {
  for (j in 1:2) if (is.na(reference[i, j])) reference[i, j] <- reference[i - 1L, j]
}
pot_rows <- !duplicated(reference[[1]])
pots <- reference[pot_rows, 1:2]
names(pots) <- c("pot_id", "modality")
if (anyNA(pots$pot_id) || anyNA(pots$modality)) stop("Incomplete reference pot metadata")
pots$pot_type <- sub(" .*", "", pots$modality)
pots$members <- strsplit(sub("^[^ ]+ +", "", pots$modality), "_", fixed = TRUE)
if (any(!pots$pot_type %in% c("mono", "mixt", "solo"))) stop("Unknown reference pot type")

# Build exactly one reference row per pot/genotype/neighbour to prevent join multiplication.
lookup <- do.call(rbind, lapply(seq_len(nrow(pots)), function(i) {
  members <- as.integer(pots$members[[i]])
  expected_n <- if (pots$pot_type[i] == "solo") 1L else 2L
  if (length(members) != expected_n || anyNA(members)) stop("Invalid reference pairing")
  do.call(rbind, lapply(unique(members), function(g) {
    neighbour <- if (length(members) == 1L) NA_integer_ else members[match(g, members) %% 2L + 1L]
    rows <- which(reference[[1]] == pots$pot_id[i] &
                    (is.na(reference[[3]]) | reference[[3]] == g))
    # Include pot-level harvest comments and only the focal genotype's comments.
    notes <- as.character(reference[[6]][rows])
    notes <- unique(notes[!is.na(notes) & nzchar(trimws(notes))])
    harvest_comment <- if (length(notes)) paste(notes, collapse = " | ") else NA_character_
    resown <- as.integer(any(reseeding_codes[rows] == "r", na.rm = TRUE))
    data.frame(pot_id = pots$pot_id[i], pot_type = pots$pot_type[i], genotype_id = g,
               neighbour_id = neighbour, resown = resown, general_comment = harvest_comment,
               check.names = FALSE)
  }))
}))
key <- function(x) paste(x$pot_id, x$pot_type, x$genotype_id, sep = "_")
if (anyDuplicated(key(lookup))) stop("Non-unique reference keys")

# Parse the diagram as displayed. Columns are ranks of pot-bearing spreadsheet
# columns across both chapelles (including offset end pots); rows are ranks of
# pot-bearing spreadsheet rows within each chapelle. These are ordinal coordinates.
layout_file <- file.path(raw, "greenhouse_layout.xlsx")
grid <- as.matrix(readxl::read_excel(layout_file, col_names = FALSE,
                                    col_types = "text", .name_repair = "minimal"))
labels <- which(grepl("^chapelle +[0-9]+$", trimws(tolower(grid))), arr.ind = FALSE)
label_cells <- arrayInd(labels, dim(grid))
if (!nrow(label_cells)) stop("No chapelle labels in layout")
label_order <- order(label_cells[, 1])
label_cells <- label_cells[label_order, , drop = FALSE]
label_rows <- label_cells[, 1]
chapelles <- as.integer(sub("^chapelle +", "", trimws(tolower(grid[label_cells]))))
numeric_cells <- which(matrix(!is.na(grid) & grepl("^[0-9]+$", trimws(grid)),
                             nrow = nrow(grid)), arr.ind = TRUE)
design <- data.frame(pot_id = as.integer(grid[numeric_cells]),
                     excel_row = numeric_cells[, 1], excel_column = numeric_cells[, 2])
block <- findInterval(design$excel_row, label_rows)
if (any(block == 0L)) stop("Pot above first chapelle label")
design$chapelle <- chapelles[block]
design$row <- as.integer(ave(design$excel_row, design$chapelle,
                            FUN = function(x) match(x, sort(unique(x)))))
design$column <- match(design$excel_column, sort(unique(design$excel_column)))
if (anyDuplicated(design$pot_id) || !setequal(design$pot_id, pots$pot_id))
  stop("Layout pots do not match identifier workbook uniquely")
if (anyDuplicated(design[c("chapelle", "row", "column")])) stop("Duplicate pot positions")
# Verify the modality next to each pot agrees with the identifier workbook.
layout_modality <- trimws(tolower(grid[cbind(design$excel_row, design$excel_column + 1L)]))
if (anyNA(layout_modality) || any(layout_modality != pots$modality[match(design$pot_id, pots$pot_id)]))
  stop("Layout and identifier workbook have different modalities")
design <- design[order(design$pot_id), ]
write_output(design, file.path(out, "greenhouse_pot_coordinates.csv"))

audit <- list()
build <- function(source, ids, traits, comments, filename) {
  # Correct only formatting and the known miuxt typo. Keep solo distinct.
  clean <- tolower(gsub("[[:space:]_]", "", ids))
  clean <- gsub("miuxt", "mixt", clean, fixed = TRUE)
  # The reference harvest note explicitly documents a label "272" instead of "372" in pot 10.
  clean[!is.na(clean) & clean == "p10mixt272"] <- "p10mixt372"
  valid <- !is.na(clean) & grepl("^p[0-9]+(mono|mixt|solo)[0-9]+$", clean)
  if (!all(valid)) stop(source, ": unrecognized sample IDs: ", paste(ids[!valid], collapse = ", "))
  parsed <- data.frame(
    pot_id = as.integer(sub("^p([0-9]+).*$", "\\1", clean)),
    pot_type = sub("^p[0-9]+(mono|mixt|solo)[0-9]+$", "\\1", clean),
    genotype_id = as.integer(sub("^p[0-9]+(mono|mixt|solo)", "", clean)))
  m <- match(key(parsed), key(lookup))
  if (anyNA(m)) stop(source, ": IDs contradict the reference: ", paste(ids[is.na(m)], collapse = ", "))
  for (n in names(traits)) {
    if (!is.numeric(traits[[n]])) stop(source, ": nonnumeric trait ", n)
    if (any(!is.na(traits[[n]]) & (!is.finite(traits[[n]]) | traits[[n]] < 0)))
      stop(source, ": invalid trait ", n)
  }
  identifiers <- lookup[m, c("pot_id", "pot_type", "genotype_id", "neighbour_id")]
  plant_id <- integer(nrow(parsed))
  for (p in unique(parsed$pot_id)) {
    ix <- which(parsed$pot_id == p)
    members <- as.integer(pots$members[[match(p, pots$pot_id)]])
    if (length(ix) > length(members)) stop(source, ": too many measurements in pot ", p)
    if (parsed$pot_type[ix[1]] == "mixt") {
      plant_id[ix] <- match(parsed$genotype_id[ix], members)
    } else plant_id[ix] <- seq_along(ix)
  }
  if (anyDuplicated(data.frame(pot_id = parsed$pot_id, plant_id = plant_id)))
    stop(source, ": duplicate plant slots within pot")
  position <- design[match(parsed$pot_id, design$pot_id), c("chapelle", "row", "column")]
  stopifnot(!anyNA(position), all(plant_id %in% 1:2))
  result <- cbind(identifiers, plant_id = plant_id, position, traits)
  result$commentaires <- comments
  result$resown <- lookup$resown[m]
  result$general_comment <- lookup$general_comment[m]
  stopifnot(nrow(result) == length(ids))
  stopifnot(all(is.na(result$neighbour_id[result$pot_type == "solo"])),
            all(result$genotype_id[result$pot_type == "mono"] == result$neighbour_id[result$pot_type == "mono"]),
            all(result$genotype_id[result$pot_type == "mixt"] != result$neighbour_id[result$pot_type == "mixt"]))
  write_output(result, file.path(processed, filename))
  saved <- read.csv2(file.path(processed, filename), fileEncoding = "UTF-8-BOM", check.names = FALSE)
  stopifnot(nrow(saved) == nrow(result), identical(names(saved), names(result)))
  stopifnot(all(saved$resown %in% c(0L, 1L)),
            identical(is.na(saved$neighbour_id), is.na(result$neighbour_id)),
            identical(saved$general_comment, result$general_comment))
  stopifnot(!any(vapply(saved, function(x) any(!is.na(x) & trimws(as.character(x)) == ""), logical(1))))
  for (n in names(traits)) stopifnot(isTRUE(all.equal(saved[[n]], result[[n]])))
  canonical <- paste0("P", parsed$pot_id, "_", parsed$pot_type, parsed$genotype_id)
  audit[[source]] <<- data.frame(source = source, source_row = seq_along(ids),
                                original_sample_id = ids, canonical_sample_id = canonical,
                                parsed, neighbour_id = lookup$neighbour_id[m],
                                plant_id = plant_id, position,
                                repeated_key = duplicated(key(parsed)) | duplicated(key(parsed), fromLast = TRUE))
  message(filename, ": ", nrow(result), " rows; all identifiers matched the reference")
}

grains <- as.data.frame(readxl::read_excel(file.path(raw, "Mosaic_pesee_grains.xlsx")))
build("Mosaic_pesee_grains.xlsx", grains[[1]],
      data.frame(poids_grains_g = grains[[2]], nb_epis = grains[[3]]), grains[[4]],
      "traits_pesee_grains.csv")
straw <- as.data.frame(readxl::read_excel(file.path(raw, "Mosaic_pesee_pailles.xlsx")))
build("Mosaic_pesee_pailles.xlsx", straw[[1]],
      data.frame(poids_pailles_g = straw[[2]]), straw[[3]], "traits_pesee_pailles.csv")
metro <- read.csv(file.path(processed, "grain_data.csv"), fileEncoding = "UTF-8", check.names = FALSE)
# grain_data.csv has no source comment field; leave commentaires empty.
build("grain_data.csv", metro$sample_id,
      data.frame(mass_totale_mesuree_g = metro[["Masse totale mesur\u00e9e"]],
                 nb_grains = metro[["Nb grains"]], PMG = metro$PMG),
      rep(NA_character_, nrow(metro)), "traits_grain_metrologie.csv")
write_output(do.call(rbind, audit), file.path(out, "trait_identifier_audit.csv"))
write_output(lookup, file.path(out, "pot_genotype_reference.csv"))
