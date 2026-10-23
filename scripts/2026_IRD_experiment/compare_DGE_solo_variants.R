# Compare competition DGE expectations with raw solo means and solo BLUEs.
# Run from project root: Rscript scripts/2026_IRD_experiment/compare_DGE_solo_variants.R
# Uses saved M3 DGE predictions; fits a separate ASReml solo fixed-genotype model.
invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
root <- dirname(dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1])))))
input <- file.path(root, 'data', '2026_IRD_experiment', 'processed_data')
out <- file.path(root, 'outputs', '2026_IRD_expermiment', 'ige_estimations')
read_csv <- function(path) read.csv2(path, fileEncoding = 'UTF-8-BOM', na.strings = 'NA')
write_csv <- function(x, filename) {
  lines <- character()
  tc <- textConnection('lines', 'w', local = TRUE)
  write.csv2(x, tc, row.names = FALSE, na = 'NA')
  close(tc)
  con <- file(file.path(out, filename), 'wb')
  on.exit(close(con))
  writeBin(as.raw(c(239,187,191)), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = '\r\n'), '\r\n'))), con)
}
selection <- read_csv(file.path(out, '01_selected_models.csv'))
dge <- read_csv(file.path(out, '03_DGE_genotype_predictions.csv'))
traits <- intersect(selection$trait[selection$selected_model == 'M3'], 'poids_grains_g')
if (!length(traits)) stop('Grain weight is not selected as M3; no comparison generated.')
d <- read_csv(file.path(input, 'traits_pesee_grains.csv'))
solo <- d[d$pot_type == 'solo' & is.finite(d$poids_grains_g), ]
stopifnot(!anyDuplicated(solo$pot_id), nrow(solo) > 0)
solo$genotype <- factor(solo$genotype_id, levels = sort(unique(solo$genotype_id)))
solo$chapelle <- factor(solo$chapelle, levels = sort(unique(solo$chapelle)))
solo$resown <- factor(solo$resown, levels = sort(unique(solo$resown)))
solo$row_group <- interaction(solo$chapelle, solo$row, drop = TRUE)
solo$column_group <- interaction(solo$chapelle, solo$column, drop = TRUE)
# BLUEs: fixed genotype, with the same estimable solo design corrections as before.
# Pot type is constant; a solo pot has no within-pot residual pair.
fit <- asreml::asreml(
  fixed = poids_grains_g ~ chapelle + resown + genotype,
  random = ~ row_group + column_group,
  residual = ~ idv(units),
  data = solo,
  na.action = asreml::na.method(y = 'fail', x = 'fail'),
  maxit = 100, workspace = '256mb'
)
if (!isTRUE(fit$converge)) fit <- update(fit, maxit = 100)
if (!isTRUE(fit$converge)) stop('Fixed-genotype solo model failed to converge.')
coefs <- summary(fit, coef = TRUE)$coef.fixed
intercept_index <- match('(Intercept)', rownames(coefs))
genotypes <- levels(solo$genotype)
genotype_indices <- match(paste0('genotype_', genotypes), rownames(coefs))
if (anyNA(c(intercept_index, genotype_indices))) stop('Cannot identify fixed coefficients.')
intercept <- coefs[intercept_index, 'solution']
blue <- coefs[genotype_indices, 'solution']
# First-level fixed effects must be zero for intercept + genotype to represent
# chapelle 1 and non-resown plants. Validate the actual ASReml parameterisation.
reference_indices <- match(c(paste0('chapelle_', levels(solo$chapelle)[1]),
                             paste0('resown_', levels(solo$resown)[1])), rownames(coefs))
stopifnot(!anyNA(reference_indices), all(coefs[reference_indices, 'solution'] == 0))
raw_means <- tapply(solo$poids_grains_g, solo$genotype, mean)
counts <- table(solo$genotype)
y <- data.frame(trait = 'poids_grains_g', genotype = genotypes,
  n_solo = as.integer(counts[genotypes]), raw_solo_mean = as.numeric(raw_means[genotypes]),
  solo_BLUE_intercept = intercept, solo_genotype_BLUE = as.numeric(blue),
  solo_BLUE_prediction = intercept + as.numeric(blue),
  reference_chapelle = levels(solo$chapelle)[1], reference_resown = levels(solo$resown)[1],
  row_column_effects = 'zero (population mean)')
x <- dge[dge$trait == 'poids_grains_g', c('trait','genotype','prediction')]
x$genotype <- as.character(x$genotype)
names(x)[3] <- 'DGE_prediction'
stopifnot(!anyDuplicated(x[c('trait','genotype')]), !anyDuplicated(y[c('trait','genotype')]))
paired <- merge(x, y, by = c('trait','genotype'), sort = FALSE)
stopifnot(nrow(paired) == length(genotypes), all(is.finite(paired$solo_BLUE_prediction)),
          all(abs(paired$solo_BLUE_prediction-paired$solo_BLUE_intercept-paired$solo_genotype_BLUE)<1e-10))
write_csv(paired, '03_DGE_solo_variants_matched_predictions.csv')
write_csv(data.frame(term = rownames(coefs), coefs, row.names = NULL), '03_solo_BLUE_fixed_effects.csv')
vc <- summary(fit)$varcomp
write_csv(data.frame(term = rownames(vc), vc, row.names = NULL), '03_solo_BLUE_variance_parameters.csv')
saveRDS(fit, file.path(out, '03_solo_BLUE_poids_grains_g_model.rds'))
results <- list()
for (variant in c('raw_mean', 'BLUE')) {
  ycol <- if (variant == 'raw_mean') 'raw_solo_mean' else 'solo_BLUE_prediction'
  good <- is.finite(paired$DGE_prediction) & is.finite(paired[[ycol]])
  xy <- paired[good, ]
  test <- cor.test(xy$DGE_prediction, xy[[ycol]], method = 'pearson')
  results[[variant]] <- data.frame(trait = 'poids_grains_g', solo_summary = variant,
    n_genotypes = nrow(xy), pearson_R = unname(test$estimate), p_nominal = test$p.value,
    note = 'Descriptive genotype-level correlation; ignores uncertainty in estimated genotype values')
  ylabel <- if (variant == 'raw_mean') 'Raw solo mean grain weight (g per plant)' else
    'Solo intercept + fixed genotype BLUE (g per plant)'
  subtitle <- sprintf('Pearson R = %.3f; nominal p = %s; n = %d genotypes',
    test$estimate, format.pval(test$p.value, digits = 3, eps = .001), nrow(xy))
  caption <- if (variant == 'raw_mean')
    'Raw solo means include all observations and receive no design correction.\nPearson p-value is exploratory and ignores uncertainty in the genotype summaries.' else
    'Solo BLUEs: chapelle 1, non-resown, random row/column effects set to zero.\nPearson p-value is exploratory and ignores uncertainty/covariance of estimated genotype values.'
  p <- ggplot(xy, aes(x = DGE_prediction, y = .data[[ycol]])) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#B95442', linewidth = .7) +
    geom_point(size = 2.7, colour = '#277DA8') +
    geom_text(aes(label = genotype), vjust = -.8, size = 3.7) +
    scale_y_continuous(expand = expansion(mult = c(.08,.15))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'),
          plot.caption = element_text(hjust = 0, size = 9)) +
    labs(title = if (variant == 'raw_mean') 'Competition DGE versus raw solo means' else
      'Competition DGE versus design-adjusted solo BLUEs', subtitle = subtitle,
      x = 'Competition-model intercept + DGE BLUP (g per plant)', y = ylabel, caption = caption)
  ggsave(file.path(out, paste0('03_DGE_vs_solo_', variant, '_poids_grains_g.png')),
         p, width = 9, height = 7, dpi = 400, bg = 'white')
}
results <- do.call(rbind, results)
write_csv(results, '03_DGE_solo_variants_correlations.csv')
print(results, row.names = FALSE)
message('Solo fixed-genotype model converged; both comparison figures generated.')
