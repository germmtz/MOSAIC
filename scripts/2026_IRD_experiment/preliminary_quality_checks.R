# Preliminary quality checks for the MOSAIC wheat competition experiment.
# Run: Rscript scripts/preliminary_quality_checks.R
# Inputs: Excel-compatible trait CSVs created by build_trait_files.R.
# Outputs: outputs/preliminary_quality_checks (400 dpi PNG + tables).
# These are exploratory, marginal comparisons, not adjusted competition models.
# Pot means avoid pseudoreplication in tests; source measurements are never changed.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
set.seed(2026)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg)) dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1])))) else getwd()
input_dir <- file.path(root, "data", "processed_data")
output_dir <- file.path(root, "outputs", "preliminary_quality_checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_trait <- function(name) read.csv2(file.path(input_dir, name),
                                      fileEncoding = "UTF-8-BOM", na.strings = "NA")
grains <- read_trait("traits_pesee_grains.csv")
straw <- read_trait("traits_pesee_pailles.csv")
metro <- read_trait("traits_grain_metrologie.csv")
traits_by_source <- list(grains = c("poids_grains_g", "nb_epis"),
                         straw = "poids_pailles_g",
                         metro = c("mass_totale_mesuree_g", "nb_grains", "PMG"))
sources <- list(grains = grains, straw = straw, metro = metro)
trait_names <- unlist(traits_by_source, use.names = FALSE)
trait_labels <- c(poids_grains_g = "Grain weight (g)", nb_epis = "Number of ears",
                  poids_pailles_g = "Straw weight (g)",
                  mass_totale_mesuree_g = "Instrument grain mass (g)",
                  nb_grains = "Number of grains", PMG = "Thousand-grain weight (g)")
mean_available <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
stars <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "***",
                           ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
p_text <- function(p) format.pval(p, digits = 2, eps = .001)
save_table <- function(x, name) write.csv2(x, file.path(output_dir, paste0(name, ".csv")),
                                          row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
save_figure <- function(p, name, width = 12, height = 8) {
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 400, bg = "white")
}
theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                  strip.background = element_rect(fill = "#EDF2F5", colour = NA),
                  strip.text = element_text(face = "bold"),
                  plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(colour = "#444444"),
                  plot.caption = element_text(hjust = 0, size = 9),
                  legend.position = "bottom"))

# Validate source keys and nonnegative numeric traits before any joins or summaries.
for (s in names(sources)) {
  d <- sources[[s]]
  stopifnot(!anyNA(d[c("pot_id", "pot_type", "genotype_id", "plant_id", "chapelle", "row", "column", "resown")]),
            !anyDuplicated(d[c("pot_id", "plant_id")]), all(d$resown %in% 0:1))
  for (t in traits_by_source[[s]]) stopifnot(is.numeric(d[[t]]),
                                           all(is.na(d[[t]]) | (is.finite(d[[t]]) & d[[t]] >= 0)))
}
long <- bind_rows(lapply(names(sources), function(s) {
  sources[[s]] %>% mutate(source = s) %>%
    pivot_longer(all_of(traits_by_source[[s]]), names_to = "trait", values_to = "value")
})) %>% mutate(trait = factor(trait, levels = trait_names))
facet_traits <- function() facet_wrap(vars(trait), scales = "free_y", ncol = 3,
                                      labeller = as_labeller(trait_labels))

# SECTION 1 -------------------------------------------------------------------
# Distributions of observed plant measurements, annotated with absolute zero counts.
# Absent rows are recorded in the coverage table.
distribution_summary <- long %>% group_by(trait) %>%
  summarise(n_rows = n(), n_observed = sum(!is.na(value)),
            missing_percent = 100 * mean(is.na(value)),
            n_zero = sum(value == 0, na.rm = TRUE),
            zero_percent = if (all(is.na(value))) NA_real_ else 100 * mean(value == 0, na.rm = TRUE),
            mean = mean_available(value), median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(annotation = sprintf("Zeros: %d\nn = %d", n_zero, n_observed))
save_table(distribution_summary, "01_distribution_summary")
design_reference <- readxl::read_excel(file.path(root, "data", "raw_data", "identifiants_echantillons_MOSAIC.xlsx"))
planned_pots <- data.frame(pot_id = design_reference[[1]],
                           pot_type = sub(" .*", "", design_reference[[2]])) %>% filter(!is.na(pot_id))
coverage <- bind_rows(lapply(names(sources), function(s) {
  planned_pots %>% left_join(count(sources[[s]], pot_id, pot_type, name = "observed_rows"),
                             by = c("pot_id", "pot_type")) %>%
    mutate(observed_rows = replace_na(observed_rows, 0L),
           source = s, expected_rows = ifelse(pot_type == "solo", 1L, 2L),
           missing_rows = expected_rows - observed_rows)
}))
save_table(coverage, "01_pot_coverage")
p_distribution <- ggplot(filter(long, !is.na(value)), aes(x = value)) +
  geom_histogram(bins = 25, fill = "#327FA5", colour = "white", linewidth = .25) +
  facet_wrap(vars(trait), scales = "free", ncol = 3, labeller = as_labeller(trait_labels)) +
  geom_label(data = distribution_summary, aes(x = Inf, y = Inf, label = annotation),
             inherit.aes = FALSE, hjust = 1.04, vjust = 1.1, size = 3.3, linewidth = 0) +
  labs(title = "Trait distributions", subtitle = "Observed plant measurements; zero values retained",
       x = "Trait value", y = "Number of measurements",
       caption = "Zeros: absolute number of zero values; n: number of nonmissing measurements. Unrecorded plants are listed in the coverage table.")
save_figure(p_distribution, "01_trait_distributions")

# Preserve the original sample identifiers and all ties for the minimum PMG.
original_metro <- read.csv(file.path(input_dir, "grain_data.csv"),
                           fileEncoding = "UTF-8", check.names = FALSE)
minimum_pmg <- min(original_metro$PMG, na.rm = TRUE)
minimum_samples <- original_metro[!is.na(original_metro$PMG) & original_metro$PMG == minimum_pmg, ]
save_table(minimum_samples, "01_minimum_PMG_samples")
message("Minimum PMG = ", minimum_pmg, "; sample(s): ", paste(minimum_samples$sample_id, collapse = ", "))

# SECTION 2 -------------------------------------------------------------------
# Marginal design associations on independent pot means. For resowing, 1 means
# at least one recorded genotype in that pot was resown. These tests do not
# isolate a causal effect or adjust for genotype, pot type, or other design terms.
pot_long <- long %>% group_by(source, trait, pot_id, pot_type, chapelle, row, column) %>%
  summarise(value = mean_available(value), resown = max(resown), .groups = "drop")
# Report counts directly from the plant-level input (before summarise overwrites value).
pot_counts <- long %>% group_by(source, trait, pot_id) %>%
  summarise(n_observed_plants = sum(!is.na(value)), .groups = "drop")
pot_long <- pot_long %>%
  left_join(pot_counts, by = c("source", "trait", "pot_id"))
save_table(pot_long, "02_pot_means")
rank_test <- function(d) {
  d <- d[is.finite(d$value) & !is.na(d$group), ]
  sizes <- table(d$group)
  if (length(sizes) < 2 || length(unique(d$value)) < 2)
    return(data.frame(n_pots = nrow(d), n_groups = length(sizes), statistic = NA_real_, p = NA_real_))
  test <- kruskal.test(value ~ group, data = d)
  data.frame(n_pots = nrow(d), n_groups = length(sizes), statistic = unname(test$statistic), p = test$p.value)
}
design_factors <- c("chapelle", "row", "column", "resown")
design_tests <- bind_rows(lapply(design_factors, function(f) {
  bind_rows(lapply(trait_names, function(t) {
    d <- filter(pot_long, trait == t) %>% mutate(group = factor(.data[[f]]))
    cbind(data.frame(factor = f, trait = t), rank_test(d))
  }))
})) %>% mutate(p_adjusted = p.adjust(p, method = "BH"), significance = stars(p_adjusted))
save_table(design_tests, "02_design_tests")
group_plot <- function(f, tests, title, caption) {
  d <- pot_long %>% filter(!is.na(value)) %>% mutate(group = factor(.data[[f]]))
  if (f == "pot_type") d$group <- factor(d$group, levels = c("solo", "mono", "mixt"))
  ann <- tests %>% mutate(trait = factor(trait, levels = trait_names),
                         label = paste0(significance, "  q = ", p_text(p_adjusted)))
  ggplot(d, aes(group, value)) +
    geom_boxplot(fill = "#9DC6D8", width = .65, outlier.shape = NA, linewidth = .4) +
    geom_point(position = position_jitter(width = .15, height = 0, seed = 2026),
               alpha = .35, size = 1.1, colour = "#254A62") + facet_traits() +
    geom_label(data = ann, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
               hjust = 1.03, vjust = 1.05, size = 3.2, linewidth = 0) +
    scale_y_continuous(expand = expansion(mult = c(.04, .18))) +
    labs(title = title, subtitle = "One point per pot; observed plant values averaged within pot",
         x = if (f == "resown") "Resown: 0 = no, 1 = at least one genotype" else f,
         y = "Pot mean trait value", caption = caption)
}
design_caption <- paste("Kruskal-Wallis omnibus tests; BH correction across 24 design tests.",
                        "\nMarginal associations, unadjusted for other factors. * q < .05; ** q < .01; *** q < .001; ns otherwise.")
for (f in design_factors) save_figure(
  group_plot(f, filter(design_tests, factor == f), paste("Design association:", f), design_caption),
  paste0("02_design_", f))

# SECTION 3 -------------------------------------------------------------------
# Compare pot types on pot means. Omnibus stars do not imply every pair differs.
# Pairwise Wilcoxon rank-sum results are exported separately with BH correction.
type_tests <- bind_rows(lapply(trait_names, function(t) {
  d <- filter(pot_long, trait == t) %>% mutate(group = factor(pot_type))
  cbind(data.frame(trait = t), rank_test(d))
})) %>% mutate(p_adjusted = p.adjust(p, method = "BH"), significance = stars(p_adjusted))
save_table(type_tests, "03_pot_type_tests")
pairwise_types <- bind_rows(lapply(trait_names, function(t) {
  bind_rows(lapply(combn(c("solo", "mono", "mixt"), 2, simplify = FALSE), function(pair) {
    a <- filter(pot_long, trait == t, pot_type == pair[1])$value
    b <- filter(pot_long, trait == t, pot_type == pair[2])$value
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    p <- if (length(a) && length(b) && length(unique(c(a, b))) > 1)
      wilcox.test(a, b, exact = FALSE)$p.value else NA_real_
    data.frame(trait = t, group_1 = pair[1], group_2 = pair[2], n_1 = length(a), n_2 = length(b), p = p)
  }))
})) %>% mutate(p_adjusted = p.adjust(p, method = "BH"), significance = stars(p_adjusted))
save_table(pairwise_types, "03_pot_type_pairwise_tests")
save_figure(group_plot("pot_type", type_tests, "Trait differences among pot types",
                       paste("Kruskal-Wallis omnibus tests; BH correction across six traits. Pairwise tests supplied separately.",
                             "\n* q < .05; ** q < .01; *** q < .001; ns otherwise. Genotype composition is not adjusted.")),
            "03_pot_type")

# SECTION 4 -------------------------------------------------------------------
# Match exact pot/type/genotype/neighbour keys, excluding mono plants: their
# within-genotype measurement order does not identify the same physical plant.
join_keys <- c("pot_id", "pot_type", "genotype_id", "neighbour_id")
a <- grains %>% filter(pot_type %in% c("solo", "mixt")) %>%
  select(all_of(join_keys), poids_grains_g)
b <- metro %>% filter(pot_type %in% c("solo", "mixt")) %>%
  select(all_of(join_keys), mass_totale_mesuree_g)
stopifnot(!anyDuplicated(a[join_keys]), !anyDuplicated(b[join_keys]))
unmatched <- bind_rows(anti_join(a, b, by = join_keys) %>% mutate(missing_from = "metrology"),
                       anti_join(b, a, by = join_keys) %>% mutate(missing_from = "grain weighing"))
save_table(unmatched, "04_unmatched_samples")
matched <- inner_join(a, b, by = join_keys) %>%
  filter(is.finite(poids_grains_g), is.finite(mass_totale_mesuree_g))
stopifnot(nrow(matched) >= 3)
agreement <- cor.test(matched$poids_grains_g, matched$mass_totale_mesuree_g, method = "pearson")
save_table(matched, "04_matched_grain_weights")
save_table(data.frame(n = nrow(matched), n_pots = n_distinct(matched$pot_id),
                      pearson_r = unname(agreement$estimate), p = agreement$p.value,
                      ci_low = agreement$conf.int[1], ci_high = agreement$conf.int[2]), "04_pearson_test")
p_agreement <- ggplot(matched, aes(poids_grains_g, mass_totale_mesuree_g)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#777777") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#254A62", fill = "#CEDFE7", linewidth = .7) +
  geom_point(aes(colour = pot_type), alpha = .65, size = 2) +
  scale_colour_manual(values = c(solo = "#D67A36", mixt = "#247CA4")) + coord_equal() +
  labs(title = "Grain weighing versus instrument mass",
       subtitle = sprintf("Pearson R = %.3f; p = %s; n = %d matched plants (%d pots)",
                          agreement$estimate, p_text(agreement$p.value), nrow(matched), n_distinct(matched$pot_id)),
       x = "Grain weight (g)", y = "Instrument grain mass (g)", colour = "Pot type",
       caption = "Solo and mixt only. Dashed line: equality. Pearson p-value and regression band are nominal:\nplants sharing a mixt pot are correlated; these are exploratory, not cluster-adjusted inference.")
save_figure(p_agreement, "04_grain_weight_agreement", 8, 7)

# SECTION 5 -------------------------------------------------------------------
# All traits are joined at POT level, not by arbitrary mono plant IDs. This includes
# mono pots safely and provides independent observational units for correlations.
# Means use available plants; unequal coverage is recorded in 02_pot_means.csv.
# Method rule: Pearson only when BOTH pot-mean marginals pass Shapiro-Wilk at .05;
# Spearman otherwise (including the discrete ear/grain counts). This is a diagnostic
# heuristic, not proof of bivariate normality. Selection is reported transparently.
wide_pots <- pot_long %>% select(pot_id, trait, value) %>%
  pivot_wider(names_from = trait, values_from = value)
normality <- bind_rows(lapply(trait_names, function(t) {
  x <- wide_pots[[t]]; x <- x[is.finite(x)]
  p <- if (length(x) >= 3 && length(x) <= 5000 && length(unique(x)) > 1) shapiro.test(x)$p.value else NA_real_
  data.frame(trait = t, n = length(x), shapiro_p = p,
             pearson_eligible = !t %in% c("nb_epis", "nb_grains") && !is.na(p) && p >= .05)
}))
save_table(normality, "05_distribution_method_checks")
correlations <- bind_rows(lapply(combn(trait_names, 2, simplify = FALSE), function(pair) {
  xy <- wide_pots[, pair]; xy <- xy[complete.cases(xy), ]
  method <- if (all(normality$pearson_eligible[match(pair, normality$trait)])) "pearson" else "spearman"
  valid <- nrow(xy) >= 3 && all(vapply(xy, function(x) length(unique(x)) > 1, logical(1)))
  test <- if (valid) cor.test(xy[[1]], xy[[2]], method = method, exact = FALSE) else NULL
  data.frame(trait_1 = pair[1], trait_2 = pair[2], n_pots = nrow(xy), method = method,
             correlation = if (valid) unname(test$estimate) else NA_real_,
             p = if (valid) test$p.value else NA_real_)
})) %>% mutate(p_adjusted = p.adjust(p, method = "BH"), significance = stars(p_adjusted))
save_table(correlations, "05_pairwise_correlations")
save_table(wide_pots, "05_matched_pot_means")
tiles <- bind_rows(correlations,
                   correlations %>% rename(trait_1 = trait_2, trait_2 = trait_1),
                   data.frame(trait_1 = trait_names, trait_2 = trait_names, correlation = 1,
                              significance = "", method = "diagonal", n_pots = NA_integer_)) %>%
  mutate(x = factor(trait_1, levels = trait_names), y = factor(trait_2, levels = rev(trait_names)),
         label = ifelse(method == "diagonal", "1.00",
                        ifelse(is.na(correlation), "NA", sprintf("%.2f%s\n%s; n=%d", correlation,
                               ifelse(significance == "ns", "", significance),
                               ifelse(method == "pearson", "P", "S"), n_pots))))
p_correlations <- ggplot(tiles, aes(x, y, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = label), size = 3.3) +
  scale_fill_gradient2(low = "#C84A4A", mid = "white", high = "#287BB5", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation", na.value = "grey90") +
  scale_x_discrete(labels = trait_labels) + scale_y_discrete(labels = trait_labels) + coord_fixed() +
  labs(title = "Pairwise trait correlations among pots", subtitle = "Available plant values averaged within pot; all pot types included",
       x = NULL, y = NULL,
       caption = "P = Pearson; S = Spearman. Pairwise complete pots; n may differ between cells.\nStars use BH-adjusted p-values across 15 pairs: * < .05; ** < .01; *** < .001. No star: q >= .05.\nPooled, exploratory associations; pot type, genotype and spatial effects are not removed.") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_figure(p_correlations, "05_trait_correlation_matrix", 11, 10)

writeLines(c("Analysis units and assumptions",
             "Section 1: observed source rows; zero values retained; no imputation or outlier removal.",
             "Sections 2-3: available plant means per pot; rank tests are marginal, not causal or adjusted.",
             "Resown at pot level: at least one genotype recorded as resown (r).",
             "BH families: 24 design tests, 6 type omnibus tests, 18 type pairwise tests, 15 correlations.",
             "Section 4: solo/mixt exact sample matches; nominal Pearson inference ignores within-pot clustering.",
             "Section 5: pot means, all pot types; incomplete pots use available plants; no mono plant matching.",
             "Missing-row coverage includes all planned pots, including entirely absent pots.",
             "Preliminary associations may reflect genotype composition, design and pot type.",
             "R reference: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html"),
           file.path(output_dir, "analysis_notes.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
message("Quality checks complete: ", output_dir)
