# IGE analysis using ASReml-R; no genetic relationship matrix.
# Run: Rscript scripts/ige_estimations.R
# Input validation only (does not require an active licence):
# Rscript scripts/ige_estimations.R --validate-only
# Outputs: outputs/ige_estimations; all CSVs use French Excel conventions.
# Confirmed specification:
# - mono + mixt only (solo has no IGE); raw plant measurements, including zeros.
# - DGE and IGE have an estimated 2x2 covariance block; genotypes are unrelated.
# - DGE:IGE is a DIRECTED interaction, NOT the covariance of DGE with IGE.
# - model comparison uses the same observations and fixed effects throughout.
# - include all available observations, including zeros, resown and commented
#   plants; resown status is adjusted as a fixed effect, not used as a filter.
# Official ASReml references:
# https://asreml.kb.vsni.co.uk/knowledge-base/lrt-asreml/
# https://cookbook.asreml.vsni.co.uk/asreml-r-examples/asr025-linear-mixed-model-correlated-residual-model-comparison-log-likelihood-ratio-test.html
# Step 4 supplies diagnostics; the request described only steps 1-3.

invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
predictions_only <- "--predictions-only" %in% args
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg)) dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1])))) else getwd()
out <- file.path(root, "outputs", "ige_estimations")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05
trait_sources <- c(poids_grains_g = "traits_pesee_grains.csv", nb_epis = "traits_pesee_grains.csv",
                   poids_pailles_g = "traits_pesee_pailles.csv", nb_grains = "traits_grain_metrologie.csv",
                   PMG = "traits_grain_metrologie.csv")
write_csv <- function(d, name) {
  lines <- character()
  tc <- textConnection("lines", "w", local = TRUE)
  write.csv2(d, tc, row.names = FALSE, na = "NA")
  close(tc)
  con <- file(file.path(out, name), "wb")
  on.exit(close(con))
  writeBin(as.raw(c(239, 187, 191)), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\r\n"), "\r\n"))), con)
}
prepare_trait <- function(trait) {
  d <- read.csv2(file.path(root, "data", "processed_data", trait_sources[[trait]]),
                  fileEncoding = "UTF-8-BOM", na.strings = "NA")
  d <- d[d$pot_type %in% c("mono", "mixt"), ]
  stopifnot(!anyNA(d[c("pot_id", "genotype_id", "neighbour_id", "plant_id", "resown")]),
            !anyDuplicated(d[c("pot_id", "plant_id")]), all(d$plant_id %in% 1:2))
  d$response <- d[[trait]]
  stopifnot(is.numeric(d$response), all(is.na(d$response) | is.finite(d$response)))
  # Complete the residual pot x plant grid. Insert missing RESPONSE slots, never
  # invented phenotypes; ASReml integrates missing responses via na.method(y='include').
  # Recover the absent genotype by swapping focal/neighbour for the observed plant.
  padded <- lapply(split(d, d$pot_id), function(p) {
    missing_slots <- setdiff(1:2, p$plant_id)
    if (length(missing_slots)) {
      stopifnot(nrow(p) == 1L, length(missing_slots) == 1L)
      added <- p[1, ]
      added$plant_id <- missing_slots
      added$genotype_id <- p$neighbour_id[1]
      added$neighbour_id <- p$genotype_id[1]
      added$response <- NA_real_
      # Read reseeding metadata from the exact other genotype, not the observed plant.
      all_source <- read.csv2(file.path(root, "data", "processed_data", "traits_pesee_grains.csv"),
                              fileEncoding = "UTF-8-BOM", na.strings = "NA")
      ref <- all_source[all_source$pot_id == added$pot_id & all_source$genotype_id == added$genotype_id, ]
      if (!nrow(ref)) stop("No metadata for missing plant in pot ", added$pot_id)
      added$resown <- ref$resown[1]
      p <- rbind(p, added)
    }
    p
  })
  d <- do.call(rbind, padded)
  d <- d[order(d$pot_id, d$plant_id), ]
  rownames(d) <- NULL
  genotypes <- sort(unique(c(d$genotype_id, d$neighbour_id)))
  stopifnot(length(genotypes) == 10L)
  d$DGE <- factor(d$genotype_id, levels = genotypes)
  d$IGE <- factor(d$neighbour_id, levels = genotypes)
  # An explicitly named factor avoids ambiguous ASReml interaction coefficient ordering.
  d$DGExIGE <- factor(paste0("g", d$genotype_id, "_n", d$neighbour_id))
  d$pot <- factor(d$pot_id)
  d$plant <- factor(d$plant_id, levels = 1:2)
  d$chapelle <- factor(d$chapelle)
  d$resown <- factor(d$resown)
  d$pot_type <- factor(d$pot_type)
  d$row_group <- interaction(d$chapelle, d$row, drop = TRUE)
  d$column_group <- interaction(d$chapelle, d$column, drop = TRUE)
  stopifnot(all(table(d$pot) == 2L), !anyDuplicated(d[c("pot", "plant")]))
  d
}
data_list <- lapply(names(trait_sources), prepare_trait)
names(data_list) <- names(trait_sources)
inventory <- do.call(rbind, lapply(names(data_list), function(t) {
  d <- data_list[[t]]
  write_csv(d, paste0("input_", t, ".csv"))
  data.frame(trait = t, pots = nlevels(d$pot), plant_slots = nrow(d),
             observed = sum(!is.na(d$response)), missing = sum(is.na(d$response)),
             zero_values = sum(d$response == 0, na.rm = TRUE), directed_pairs = nlevels(d$DGExIGE))
}))
write_csv(inventory, "input_inventory.csv")

# Step 3 can be refreshed independently using saved selected models:
# Rscript scripts/ige_estimations.R --predictions-only
export_predictions <- function(selection, fitted_models = NULL) {
  selected_traits <- selection$trait[selection$selected_model == "M3"]
  solo_predictions <- data.frame()
  dge_predictions <- data.frame()
  mono_fixed_effects <- setNames(numeric(), character())
  predictions <- data.frame(trait = character(), genotype = character(), neighbour = character(),
                             interaction = numeric(), interaction_type = character())
  for (t in selected_traits) {
    m <- selection$selected_model[selection$trait == t]
    fit <- if (is.null(fitted_models)) readRDS(file.path(out, paste0(t, "_", m, ".rds"))) else fitted_models[[t]][[m]]
    blup <- function(term, level) {
      b <- fit$coefficients$random
      labels <- if (!is.null(rownames(b))) rownames(b) else names(b)
      index <- match(paste0(term, "_", level), labels)
      if (anyNA(index)) stop("Cannot identify saved BLUPs for ", term)
      as.numeric(b[index])
    }
    genotypes <- levels(data_list[[t]]$DGE)
    competition_fixed <- fit$coefficients$fixed
    intercept_index <- match("(Intercept)", rownames(competition_fixed))
    if (is.na(intercept_index)) stop("Cannot identify competition intercept")
    competition_intercept <- as.numeric(competition_fixed[intercept_index])
    mono_index <- match("pot_type_mono", rownames(competition_fixed))
    mixt_index <- match("pot_type_mixt", rownames(competition_fixed))
    if (anyNA(c(mono_index, mixt_index))) stop("Cannot identify mono-versus-mixt fixed contrast")
    mono_fixed_effects[t] <- as.numeric(competition_fixed[mono_index] - competition_fixed[mixt_index])
    dge_blup <- blup("DGE", genotypes)
    dge_predictions <- rbind(dge_predictions, data.frame(trait = t, genotype = genotypes,
      intercept = competition_intercept, DGE_BLUP = dge_blup,
      prediction = competition_intercept + dge_blup,
      prediction_formula = 'intercept + DGE',
      reference_chapelle = levels(data_list[[t]]$chapelle)[1],
      reference_pot_type = levels(data_list[[t]]$pot_type)[1],
      reference_resown = levels(data_list[[t]]$resown)[1],
      other_random_effects = "zero (population mean)"))
    pairs <- expand.grid(genotype = genotypes, neighbour = genotypes, stringsAsFactors = FALSE)
    effect <- blup("IGE", pairs$neighbour)
    if (m == "M3") effect <- effect + blup("DGExIGE", paste0("g", pairs$genotype, "_n", pairs$neighbour))
    predictions <- rbind(predictions, data.frame(trait = t, genotype = pairs$genotype,
      neighbour = pairs$neighbour, interaction = effect,
      interaction_type = if (m == "M2") "IGE only" else "IGE + DGE:IGE"))

    # Solo performance is estimated from SOLO observations, not from the IGE model.
    # Same applicable design corrections: chapelle/resown fixed, rows and columns
    # nested within chapelle random. Pot type is constant; no within-pot residual
    # correlation is estimable for single plants, so residuals are independent.
    solo <- read.csv2(file.path(root, "data", "processed_data", trait_sources[[t]]),
                       fileEncoding = "UTF-8-BOM", na.strings = "NA")
    solo <- solo[solo$pot_type == "solo" & !is.na(solo[[t]]), ]
    stopifnot(!anyDuplicated(solo$pot_id), nrow(solo) > 0)
    solo$response <- solo[[t]]
    solo$DGE <- factor(solo$genotype_id, levels = genotypes)
    solo$chapelle <- factor(solo$chapelle, levels = sort(unique(solo$chapelle)))
    solo$resown <- factor(solo$resown, levels = sort(unique(solo$resown)))
    solo$row_group <- interaction(solo$chapelle, solo$row, drop = TRUE)
    solo$column_group <- interaction(solo$chapelle, solo$column, drop = TRUE)
    solo_model_path <- file.path(out, paste0("03_solo_", t, "_model.rds"))
    if (predictions_only && file.exists(solo_model_path)) {
      solo_fit <- readRDS(solo_model_path)
    } else {
      solo_fit <- asreml::asreml(fixed = response ~ chapelle + resown,
      random = ~ row_group + column_group + DGE, residual = ~idv(units),
      data = solo, na.action = asreml::na.method(y = "fail", x = "fail"),
      maxit = 100, workspace = "256mb")
    }
    if (!isTRUE(solo_fit$converge)) solo_fit <- update(solo_fit, maxit = 100)
    if (!isTRUE(solo_fit$converge)) stop("Solo model did not converge for ", t)
    fixed_coef <- solo_fit$coefficients$fixed
    intercept_index <- match("(Intercept)", rownames(fixed_coef))
    if (is.na(intercept_index)) stop("Cannot identify solo intercept")
    intercept <- as.numeric(fixed_coef[intercept_index])
    random_coef <- solo_fit$coefficients$random
    genotype_index <- match(paste0("DGE_", genotypes), rownames(random_coef))
    if (anyNA(genotype_index)) stop("Cannot identify solo genotype BLUPs")
    solo_blup <- as.numeric(random_coef[genotype_index])
    solo_predictions <- rbind(solo_predictions, data.frame(trait = t, genotype = genotypes,
      intercept = intercept, genotype_BLUP = solo_blup, prediction = intercept + solo_blup,
      reference_chapelle = levels(solo$chapelle)[1], reference_resown = levels(solo$resown)[1],
      row_column_effects = "zero (population mean)",
      n_observed = as.integer(table(solo$DGE)[genotypes])))
    saveRDS(solo_fit, file.path(out, paste0("03_solo_", t, "_model.rds")))
    write_csv(cbind(data.frame(trait = t, term = rownames(summary(solo_fit)$varcomp)),
                    summary(solo_fit)$varcomp), paste0("03_solo_", t, "_variance_parameters.csv"))
  }
  # Delete only step-3 artifacts, after checking all saved fits can be extracted.
  prediction_dir <- file.path(out, "genotype_predictions")
  dir.create(prediction_dir, showWarnings = FALSE)
  old_files <- c(file.path(out, c("03_genotype_neighbour_interactions.csv", "03_same_genotype_interactions.csv", "03_solo_genotype_predictions.csv")),
                 list.files(prediction_dir, pattern = "\\.csv$", full.names = TRUE))
  old_files <- old_files[file.exists(old_files)]
  checked_root <- paste0(normalizePath(out, winslash = "/"), "/")
  stopifnot(all(startsWith(normalizePath(old_files, winslash = "/"), checked_root)))
  if (length(old_files) && any(!file.remove(old_files))) stop("Could not remove old prediction files")
  between <- predictions[predictions$genotype != predictions$neighbour, ]
  between$genetic_interaction <- between$interaction
  between$mono_fixed_effect <- rep(0, nrow(between))
  same <- predictions[predictions$genotype == predictions$neighbour, ]
  same$genetic_interaction <- same$interaction
  same$mono_fixed_effect <- unname(mono_fixed_effects[same$trait])
  same$interaction <- same$genetic_interaction + same$mono_fixed_effect
  same$interaction_type <- rep("mono fixed effect + IGE + DGE:IGE", nrow(same))
  # This is an adjusted interaction contribution, not a full phenotype:
  # the intercept and focal DGE remain excluded.
  # One table containing all 100 directed pairs per selected trait.
  # Mono rows include their fixed pot-type contrast; mixt rows have adjustment zero.
  all_pairs <- rbind(between, same)
  all_pairs <- all_pairs[order(all_pairs$trait, as.integer(all_pairs$genotype),
                              as.integer(all_pairs$neighbour)), ]
  stopifnot(nrow(all_pairs) == 100L * length(selected_traits),
            !anyDuplicated(all_pairs[c("trait", "genotype", "neighbour")]))
  all_pairs$prediction_formula <- all_pairs$interaction_type
  all_pairs <- all_pairs[c('trait', 'genotype', 'neighbour', 'interaction', 'prediction_formula', 'interaction_type')]
  write_csv(all_pairs, "03_all_genotype_interactions.csv")
  if (dir.exists(prediction_dir) && !length(list.files(prediction_dir, all.files = TRUE, no.. = TRUE)))
    unlink(prediction_dir, recursive = TRUE)
  if (!nrow(solo_predictions)) solo_predictions <- data.frame(trait = character(), genotype = character(),
    intercept = numeric(), genotype_BLUP = numeric(), prediction = numeric(),
    reference_chapelle = character(), reference_resown = character(), row_column_effects = character(), n_observed = integer())
  write_csv(solo_predictions, "03_solo_genotype_predictions.csv")
  if (!nrow(dge_predictions)) dge_predictions <- data.frame(trait = character(), genotype = character(),
    intercept = numeric(), DGE_BLUP = numeric(), prediction = numeric(), reference_chapelle = character(),
    reference_pot_type = character(), reference_resown = character(), other_random_effects = character())
  write_csv(dge_predictions, "03_DGE_genotype_predictions.csv")

  # Match genotype estimates, not plant observations. Pearson inference is nominal:
  # BLUPs are estimated/shrunken quantities with shared parameter uncertainty.
  paired <- merge(dge_predictions[c("trait", "genotype", "prediction")],
                  solo_predictions[c("trait", "genotype", "prediction")],
                  by = c("trait", "genotype"), suffixes = c("_DGE", "_solo"))
  write_csv(paired, "03_DGE_solo_matched_predictions.csv")
  correlation_results <- list()
  for (t in selected_traits) {
    xy <- paired[paired$trait == t & complete.cases(paired), ]
    valid <- nrow(xy) >= 3 && sd(xy$prediction_DGE) > 0 && sd(xy$prediction_solo) > 0
    test <- if (valid) cor.test(xy$prediction_DGE, xy$prediction_solo, method = "pearson") else NULL
    correlation_results[[t]] <- data.frame(trait = t, n_genotypes = nrow(xy),
      pearson_R = if (valid) unname(test$estimate) else NA_real_,
      p_nominal = if (valid) test$p.value else NA_real_,
      solo_prediction_range = diff(range(xy$prediction_solo)),
      note = "Descriptive correlation of BLUP-based predictions; ignores estimation uncertainty")
    if (nrow(xy)) {
      subtitle <- if (valid) sprintf("Pearson R = %.3f; nominal p = %s; n = %d genotypes",
        test$estimate, format.pval(test$p.value, digits = 3, eps = .001), nrow(xy)) else "Correlation not estimable"
      p <- ggplot2::ggplot(xy, ggplot2::aes(prediction_DGE, prediction_solo)) +
        ggplot2::geom_point(size = 2.5, colour = "#277DA8") +
        ggplot2::geom_text(ggplot2::aes(label = genotype), vjust = -0.7, size = 3.5, check_overlap = TRUE) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#B95442", linewidth = .6) +
        ggplot2::scale_y_continuous(labels = function(x) sprintf("%.6f", x), expand = ggplot2::expansion(mult = .15)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = paste("DGE versus solo genotype expectations:", t), subtitle = subtitle,
          x = "Competition-model intercept + DGE BLUP (trait units)",
          y = "Solo-model intercept + genotype BLUP (trait units)",
          caption = paste0("Reference fixed-effect levels are recorded in the prediction tables. Pearson p-value is exploratory.",
            if (diff(range(xy$prediction_solo)) < 1e-4 * max(1, abs(mean(xy$prediction_solo))))
              "\nSolo predictions are nearly constant: the displayed differences are magnified and strongly shrunken." else ""))
      ggplot2::ggsave(file.path(out, paste0("03_DGE_vs_solo_", t, ".png")), p, width = 9, height = 7, dpi = 400)
    }
  }
  correlation_table <- if (length(correlation_results)) do.call(rbind, correlation_results) else
    data.frame(trait = character(), n_genotypes = integer(), pearson_R = numeric(), p_nominal = numeric(),
               solo_prediction_range = numeric(), note = character())
  write_csv(correlation_table, "03_DGE_solo_correlations.csv")
  stopifnot(nrow(between) == 90L * length(selected_traits), !anyNA(between$interaction))
  notes_path <- file.path(out, "analysis_notes.txt")
  if (file.exists(notes_path)) {
    notes <- readLines(notes_path, warn = FALSE)
    notes <- notes[!grepl("^M1 predictions|^Step 3 predictions", notes)]
    notes <- notes[!grepl("^One CSV per focal|^Solo predictions|^DGE predictions|^DGE-solo correlation", notes)]
    writeLines(c(notes, "Step 3 predictions include only traits selected as M3; grouped interaction, solo and DGE tables; no per-genotype files.",
      "Solo predictions: solo-only ASReml; intercept + genotype BLUP at stated reference chapelle/resown; random row/column effects zero.",
      "DGE predictions: competition-model intercept + DGE BLUP; reference chapelle, pot type and resown levels recorded in CSV.",
      "DGE-solo correlation is descriptive Pearson across matched genotypes; nominal p ignores BLUP estimation uncertainty."), notes_path)
  }
  message("Step 3 exported ", length(selected_traits), " trait(s), ", nrow(between), " between-genotype predictions.")
}

run_analysis <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) stop("ASReml-R is not installed")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  # Same genotype levels/order for both effects: unstructured 2x2 covariance
  # across DGE/IGE, identity across genotypes (no GRM).
  random_formulas <- list(M1 = ~ row_group + column_group + DGE,
                          M2 = ~ row_group + column_group + str(~DGE + IGE, ~corgh(2):id(DGE)),
                          M3 = ~ row_group + column_group + str(~DGE + IGE, ~corgh(2):id(DGE)) + DGExIGE)
  fit_model <- function(d, model) {
    # corv = one marginal residual variance + a freely estimated correlation;
    # distinct pots are independent. Do not also add a random pot intercept:
    # that would confound the shared-pot covariance with this residual structure.
    fit <- asreml::asreml(fixed = response ~ chapelle + pot_type + resown,
                          random = random_formulas[[model]],
                          residual = ~ id(pot):corv(plant),
                          data = d, na.action = asreml::na.method(y = "include", x = "fail"),
                          maxit = 100, workspace = "256mb")
    if (!isTRUE(fit$converge)) fit <- update(fit, maxit = 100)
    if (!isTRUE(fit$converge)) stop(model, " failed to converge")
    # Correlated effects may converge to a poor boundary solution. Restart from
    # negative, zero and positive correlations with nonzero genetic variances;
    # retain the converged solution having the highest REML likelihood.
    start_log <- data.frame(start = "default", converged = fit$converge, loglik = fit$loglik)
    if (model != "M1") {
      base_fit <- fit
      for (rho in c(-.8, 0, .8)) {
        gp <- base_fit$G.param
        for (term in names(gp)) for (j in seq_along(gp[[term]])) {
          block <- gp[[term]][[j]]
          if (!is.null(block$con)) block$con[block$con == "B"] <- "P"
          gp[[term]][[j]] <- block
        }
        gp[["DGE+IGE"]][["2"]]$initial[] <- c(rho, .15 * var(d$response, na.rm = TRUE), .05 * var(d$response, na.rm = TRUE))
        gp[["DGE+IGE"]][["2"]]$con <- c("U", "P", "P")
        candidate <- try(asreml::asreml(fixed = response ~ chapelle + pot_type + resown,
          random = random_formulas[[model]], residual = ~id(pot):corv(plant),
          data = d, G.param = gp, R.param = base_fit$R.param,
          na.action = asreml::na.method(y = "include", x = "fail"),
          maxit = 100, workspace = "256mb", trace = FALSE), silent = TRUE)
        valid <- !inherits(candidate, "try-error") && isTRUE(candidate$converge)
        start_log <- rbind(start_log, data.frame(start = paste0("rho=", rho), converged = valid,
          loglik = if (valid) candidate$loglik else NA_real_))
        if (valid && candidate$loglik > fit$loglik) fit <- candidate
      }
    }
    attr(fit, "starting_value_checks") <- start_log
    vc <- summary(fit)$varcomp
    if ("bound" %in% names(vc) && any(as.character(vc$bound) %in% c("S", "C")))
      stop(model, ": singular or constrained variance parameters; inspect model before LRT")
    fit
  }
  # STEP 1: IGE detection ------------------------------------------------------
  # REML LRTs are valid here because fixed effects and response rows are identical.
  # M1 -> M2 adds IGE variance AND DGE-IGE covariance: approximate null
  # 0.5 chi-square(1) + 0.5 chi-square(2), not the one-variance mixture.
  # M2 -> M3 adds one variance: 0.5 chi-square(0) + 0.5 chi-square(1).
  # This is an approximate variance-component test, particularly with ten genotypes.
  # Raw p-values drive prespecified hierarchical selection; BH supplied separately.
  fitted_models <- list(); lrt_rows <- list(); selections <- list(); vc_rows <- list()
  compare_models <- function(reduced, full, trait, reduced_name, full_name) {
    delta <- full$loglik - reduced$loglik
    if (!is.finite(delta) || delta < -1e-4) stop("Invalid nested likelihood ordering: ", trait)
    lr <- max(0, 2 * delta)
    joint <- reduced_name == "M1"
    p <- if (lr == 0) 1 else if (joint)
      .5 * pchisq(lr, 1, lower.tail = FALSE) + .5 * pchisq(lr, 2, lower.tail = FALSE) else
        .5 * pchisq(lr, 1, lower.tail = FALSE)
    data.frame(trait = trait, reduced_model = reduced_name, full_model = full_name,
               loglik_reduced = reduced$loglik, loglik_full = full$loglik,
               LR = lr, added_parameters = if (joint) 2L else 1L, p = p,
               reference = if (joint) "0.5 chi-square(1) + 0.5 chi-square(2)" else "0.5 chi-square(0) + 0.5 chi-square(1)")
  }
  for (t in names(data_list)) {
    models <- lapply(names(random_formulas), function(m) {
      message("Fitting ", t, " / ", m)
      fit <- fit_model(data_list[[t]], m)
      saveRDS(fit, file.path(out, paste0(t, "_", m, ".rds")))
      vc <- summary(fit)$varcomp
      vc_rows[[paste(t, m)]] <<- cbind(data.frame(trait = t, model = m, term = rownames(vc)), vc)
      fit
    })
    names(models) <- names(random_formulas)
    fitted_models[[t]] <- models
    l12 <- compare_models(models$M1, models$M2, t, "M1", "M2")
    l23 <- compare_models(models$M2, models$M3, t, "M2", "M3")
    best <- if (l23$p < alpha) "M3" else if (l12$p < alpha) "M2" else "M1"
    selections[[t]] <- data.frame(trait = t, selected_model = best, alpha = alpha,
                                  rule = "M3 if M3>M2; otherwise M2 if M2>M1; otherwise M1")
    lrt_rows[[t]] <- rbind(l12, l23)
  }
  lrt_table <- do.call(rbind, lrt_rows)
  lrt_table$p_BH <- p.adjust(lrt_table$p, method = "BH")
  selection <- do.call(rbind, selections)
  lrt_table$selected_model <- selection$selected_model[match(lrt_table$trait, selection$trait)]
  write_csv(lrt_table, "01_likelihood_ratio_tests.csv")
  write_csv(selection, "01_selected_models.csv")
  write_csv(do.call(rbind, vc_rows), "01_all_variance_parameters.csv")

  # STEP 2: variance decomposition ---------------------------------------------
  # Marginal conditional phenotypic variance = genetic + row + column + residual
  # variances, conditional on fixed effects. Correlation rho is not a variance:
  # it changes off-diagonal residual covariance, not the individual denominator.
  # DGE/IGE covariance contributes 2*Cov(D,I) only when focal and neighbour are
  # the same genotype. Average over the observed mono/mixt composition, retaining
  # a signed covariance contribution (do not disguise it as positive variance).
  variance_rows <- list(); covariance_rows <- list()
  component_value <- function(vc, term, required = TRUE) {
    idx <- which(rownames(vc) == term)
    if (length(idx) == 0L && !required) return(0)
    if (length(idx) != 1L) stop("Cannot identify variance component: ", term)
    as.numeric(vc[idx, "component"])
  }
  genetic_block <- function(vc, model) {
    if (model == "M1") return(c(DGE = component_value(vc, "DGE"), IGE = 0, covariance = 0, correlation = NA_real_))
    vd <- component_value(vc, "DGE+IGE!corgh(2)_1")
    vi <- component_value(vc, "DGE+IGE!corgh(2)_2")
    ci <- grep("^DGE\\+IGE!.*\\.cor$", rownames(vc))
    if (length(ci) != 1L) stop("Cannot identify DGE-IGE correlation")
    rho <- as.numeric(vc[ci, "component"])
    c(DGE = vd, IGE = vi, covariance = rho * sqrt(vd * vi), correlation = rho)
  }
  for (t in names(data_list)) for (m in names(random_formulas)) {
    g <- genetic_block(summary(fitted_models[[t]][[m]])$varcomp, m)
    covariance_rows[[paste(t, m)]] <- data.frame(trait = t, model = m,
      DGE_variance = g[["DGE"]], IGE_variance = g[["IGE"]],
      DGE_IGE_covariance = g[["covariance"]], DGE_IGE_correlation = g[["correlation"]],
      near_correlation_boundary = !is.na(g[["correlation"]]) && abs(g[["correlation"]]) > .98)
  }
  write_csv(do.call(rbind, covariance_rows), "02_DGE_IGE_covariances.csv")
  for (t in names(data_list)) {
    selected <- selection$selected_model[selection$trait == t]
    fit <- fitted_models[[t]][[selected]]
    vc <- summary(fit)$varcomp
    # corv variance naming is engine-dependent; require exactly one residual VAR row.
    residual_index <- grep("!var$", rownames(vc), ignore.case = TRUE)
    if (length(residual_index) != 1L) stop("Cannot identify corv residual variance; inspect ", t, " fit")
    g <- genetic_block(vc, selected)
    d <- data_list[[t]]
    mono_fraction <- mean(d$genotype_id[!is.na(d$response)] == d$neighbour_id[!is.na(d$response)])
    components <- c(DGE = g[["DGE"]], IGE = g[["IGE"]],
                    `DGE-IGE covariance contribution` = 2 * mono_fraction * g[["covariance"]],
                    `DGE:IGE` = component_value(vc, "DGExIGE", selected == "M3"),
                    Row = component_value(vc, "row_group"),
                    Column = component_value(vc, "column_group"),
                    Residual = as.numeric(vc[residual_index, "component"]))
    if (any(!is.finite(components)) || sum(components) <= 0 ||
        any(components[names(components) != "DGE-IGE covariance contribution"] < 0)) stop("Invalid variance estimates for ", t)
    variance_rows[[t]] <- data.frame(trait = t, model = selected, component = names(components),
                                    variance = as.numeric(components),
                                    percent = 100 * components / sum(components),
                                    observed_mono_fraction = mono_fraction)

  }
  variance_table <- do.call(rbind, variance_rows)
  write_csv(variance_table, "02_variance_decomposition.csv")
  variance_table$component <- factor(variance_table$component,
                                     levels = c("DGE", "IGE", "DGE:IGE", "DGE-IGE covariance contribution", "Row", "Column", "Residual"))
  p <- ggplot2::ggplot(variance_table, ggplot2::aes(trait, percent, fill = component)) +
    ggplot2::geom_col(width = .7) +
    ggplot2::scale_fill_manual(values = c("#277DA8", "#E69F00", "#8A62AA", "#B95442", "#B9C7C9", "#D4DDDE", "#EEEEEE")) +
    ggplot2::labs(x = NULL, y = "Marginal phenotypic variance (%)", fill = "Component",
                  caption = "Conditional on fixed effects; signed contributions sum to 100%.\nCovariance contribution = 2 x observed mono fraction x Cov(DGE, IGE); negative values plot below zero.") +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  ggplot2::ggsave(file.path(out, "02_variance_decomposition.png"), p, width = 10, height = 6, dpi = 400)
  export_predictions(selection, fitted_models)
  # STEP 4: diagnostics and reproducibility -------------------------------------
  diagnostics <- do.call(rbind, lapply(names(data_list), function(t) {
    do.call(rbind, lapply(names(random_formulas), function(m) {
      f <- fitted_models[[t]][[m]]
      data.frame(trait = t, model = m, converged = isTRUE(f$converge), loglik = f$loglik,
                  variance_parameters = length(f$vparameters))
    }))
  }))
  write_csv(diagnostics, "04_fit_diagnostics.csv")
  starts <- do.call(rbind, lapply(names(fitted_models), function(t)
    do.call(rbind, lapply(names(fitted_models[[t]]), function(m)
      cbind(data.frame(trait = t, model = m), attr(fitted_models[[t]][[m]], "starting_value_checks"))))))
  write_csv(starts, "04_starting_value_checks.csv")
  residual_data <- do.call(rbind, lapply(names(data_list), function(t) {
    m <- selection$selected_model[selection$trait == t]
    f <- fitted_models[[t]][[m]]
    d <- data_list[[t]]
    e <- as.numeric(f$residuals)
    fitted <- as.numeric(f$linear.predictors)
    stopifnot(length(e) == nrow(d), length(fitted) == nrow(d))
    observed <- !is.na(d$response)
    stopifnot(max(abs(d$response[observed] - fitted[observed] - e[observed])) < 1e-6)
    data.frame(trait = t, model = m, pot_id = d$pot_id[observed], plant_id = d$plant_id[observed],
                fitted = fitted[observed], residual = e[observed])
  }))
  write_csv(residual_data, "04_selected_model_residuals.csv")
  qq <- ggplot2::ggplot(residual_data, ggplot2::aes(sample = residual)) +
    ggplot2::stat_qq(alpha = .5, size = 1) + ggplot2::stat_qq_line(colour = "#B95442") +
    ggplot2::facet_wrap(~trait, scales = "free", ncol = 3) + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "Theoretical normal quantiles", y = "Conditional residual", title = "Selected-model residual Q-Q plots")
  ggplot2::ggsave(file.path(out, "04_residual_qq.png"), qq, width = 11, height = 7, dpi = 300)
  rvf <- ggplot2::ggplot(residual_data, ggplot2::aes(fitted, residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "#B95442") +
    ggplot2::geom_point(alpha = .4, size = 1) + ggplot2::facet_wrap(~trait, scales = "free", ncol = 3) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::labs(x = "Fitted value", y = "Conditional residual")
  ggplot2::ggsave(file.path(out, "04_residual_vs_fitted.png"), rvf, width = 11, height = 7, dpi = 300)
  writeLines(c("Model specification and interpretation",
               "All available observations retained: no comment-based, outlier or resown exclusions; zero values retained.",
               "Resown status enters as a fixed effect in competition and solo models; mono flags remain pot-genotype-level annotations.",
               "Mono and mixt only. DGE/IGE have an estimated 2x2 covariance block, identity across genotypes; no GRM.",
               "Fixed: chapelle + pot_type + resown. Random nuisance: row and column groups within chapelle.",
               "Residual: id(pot):corv(plant), one marginal variance and within-pot residual correlation.",
               "DGE:IGE is a directed pair interaction, not DGE-IGE covariance.",
               "M1-M2 jointly tests IGE variance and covariance: 50:50 chi-square(1)/chi-square(2) approximation.",
               "M2-M3 tests the additional interaction variance: 50:50 chi-square(0)/chi-square(1) approximation.",
               "With ten genotypes and nuisance components near boundaries, asymptotic p-values warrant caution.",
               "Selection uses raw p<0.05, with interaction tested first; supplementary BH covers ten tests.",
               "Variance percentages include signed 2*mono_fraction*Cov(DGE,IGE), plus genetic/environmental/residual variances.",
               "Across different unrelated genotypes DGE-IGE covariance has zero diagonal contribution; across mono plants it contributes.",
               "Near-boundary DGE-IGE correlations are flagged in 02_DGE_IGE_covariances.csv; asymptotic inference may be fragile.",
               "Prediction genotype = focal recipient, neighbour = actor producing the IGE; values are BLUP deviations.",
               "Step 3 predictions include only traits selected as M3; grouped interaction, solo and DGE tables; no per-genotype files.",
               "DGE predictions: competition-model intercept + DGE BLUP; reference fixed-effect levels recorded in CSV.",
               "DGE-solo correlation is descriptive Pearson across matched genotypes; nominal p ignores BLUP estimation uncertainty.",
               "Solo predictions: solo-only ASReml; intercept + genotype BLUP at stated reference chapelle/resown; random row/column effects zero.",
               "Gaussian models were requested for all traits; inspect residual plots, particularly count traits.",
               "No licence codes are stored in analysis scripts or outputs."),
             file.path(out, "analysis_notes.txt"))
  capture.output(sessionInfo(), file = file.path(out, "04_session_info.txt"))
}

if (predictions_only) {
  selection <- read.csv2(file.path(out, "01_selected_models.csv"), fileEncoding = "UTF-8-BOM")
  export_predictions(selection)
  writeLines("Step 3 refreshed using saved competition fits; solo fits reused when available or fitted separately.",
             file.path(out, "analysis_status.txt"))
} else if (validate_only) {
  writeLines(c("Input validation succeeded. Model fits were NOT attempted.",
               "Model specification: mono/mixt, correlated DGE/IGE, directed interaction.",
               "An active ASReml-R licence is required to run the full analysis."),
             file.path(out, "validation_status.txt"))
  message("Input validation complete; no model fits attempted.")
} else {
  tryCatch(run_analysis(), error = function(e) {
    writeLines(c("Analysis did not complete. No uncomputed estimates should be interpreted.",
                 conditionMessage(e)), file.path(out, "analysis_status.txt"))
    stop(e)
  })
  writeLines("Analysis completed successfully.", file.path(out, "analysis_status.txt"))
}
