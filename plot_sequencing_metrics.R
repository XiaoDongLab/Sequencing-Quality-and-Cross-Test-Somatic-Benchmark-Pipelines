suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
  library(ggbeeswarm)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- ifelse(length(args) >= 1, args[1], "./bridge_tsv")
out_dir <- ifelse(length(args) >= 2, args[2], "./plots_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

depth_curve_path <- file.path(in_dir, "depth_curve.tsv")
lorenz_path      <- file.path(in_dir, "lorenz_curve.tsv")
dist_path        <- file.path(in_dir, "depth_distribution.tsv")
chrom_path       <- file.path(in_dir, "chrom_coverage_1Mb.tsv")
vaf_path         <- file.path(in_dir, "vaf_values.tsv")

boot_ci_mean <- function(x, B = 500, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  m <- replicate(B, mean(sample(x, replace = TRUE)))
  alpha <- (1 - conf) / 2
  quantile(m, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
}

calc_gini_from_lorenz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) == 0) return(NA_real_)
  o <- order(x, y)
  x <- x[o]; y <- y[o]
  if (x[1] > 0 || y[1] > 0) { x <- c(0, x); y <- c(0, y) }
  if (x[length(x)] < 1 || y[length(y)] < 1) { x <- c(x, 1); y <- c(y, 1) }
  dx <- diff(x)
  auc <- sum(dx * (y[-1] + y[-length(y)]) / 2)
  1 - 2 * auc
}

theme_base <- theme_classic() +
  theme(
    axis.title = element_text(color = "black"),
    axis.text  = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.title = element_blank()
  )

depth_curve <- read_tsv(depth_curve_path, show_col_types = FALSE) %>%
  mutate(
    reads_number = as.numeric(reads_number),
    coverage_frac = as.numeric(coverage_frac)
  ) %>%
  filter(is.finite(reads_number), is.finite(coverage_frac))

depth_curve_sum <- depth_curve %>%
  group_by(reads_number) %>%
  summarise(
    n = n_distinct(sample_id),
    mean_cov = mean(coverage_frac, na.rm = TRUE),
    ci = list(boot_ci_mean(coverage_frac)),
    .groups = "drop"
  ) %>%
  mutate(
    ci_low = vapply(ci, `[[`, numeric(1), 1),
    ci_high = vapply(ci, `[[`, numeric(1), 2)
  )

p1 <- ggplot(depth_curve_sum, aes(x = reads_number, y = mean_cov)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.7) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "Reads-versus-coverage curve (mean with bootstrap 95% CI)",
    x = "Downsampled reads",
    y = "Genome coverage"
  ) +
  theme_base

ggsave(file.path(out_dir, "reads_coverage_curve.pdf"), p1, width = 6.5, height = 4.2)

a <- read_tsv(lorenz_path, show_col_types = FALSE) %>%
  mutate(
    x_frac = as.numeric(x_frac),
    y_frac = as.numeric(y_frac)
  ) %>%
  filter(is.finite(x_frac), is.finite(y_frac)) %>%
  group_by(sample_id) %>%
  arrange(x_frac, .by_group = TRUE) %>%
  ungroup()

lorenz_binned <- a %>%
  mutate(x_bin = pmin(199L, pmax(0L, as.integer(floor(x_frac * 200))))) %>%
  group_by(x_bin) %>%
  summarise(
    x = mean(x_frac),
    y_med = median(y_frac),
    y_lo = quantile(y_frac, 0.1),
    y_hi = quantile(y_frac, 0.9),
    .groups = "drop"
  )

gini_tbl <- a %>%
  group_by(sample_id) %>%
  summarise(gini = calc_gini_from_lorenz(x_frac, y_frac), .groups = "drop")

p2a <- ggplot(lorenz_binned, aes(x = x, y = y_med)) +
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi), alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "Lorenz curve (median with 10-90% band across samples)",
    x = "Sequenced bases fraction",
    y = "Genome covered fraction"
  ) +
  theme_base

p2b <- ggplot(gini_tbl, aes(x = "", y = gini)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.35) +
  stat_summary(fun = median, geom = "point", size = 2.2) +
  labs(title = "Gini coefficient across samples", x = NULL, y = "Gini (from Lorenz)") +
  theme_base +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- p2a / p2b + plot_layout(heights = c(2.2, 1.4))
ggsave(file.path(out_dir, "lorenz_curve_and_gini.pdf"), p2, width = 6.5, height = 7.2)
write_tsv(gini_tbl, file.path(out_dir, "gini_by_sample.tsv"))

dist <- read_tsv(dist_path, show_col_types = FALSE) %>%
  mutate(depth = as.integer(depth), count = as.numeric(count)) %>%
  filter(is.finite(depth), is.finite(count), depth >= 0, depth <= 100)

grid <- tidyr::expand_grid(sample_id = unique(dist$sample_id), depth = 0:100)

dist2 <- grid %>%
  left_join(dist, by = c("sample_id", "depth")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(sample_id) %>%
  mutate(total = sum(count), count_frac = ifelse(total > 0, count / total, NA_real_)) %>%
  ungroup()

dist_sum <- dist2 %>%
  group_by(depth) %>%
  summarise(
    y_med = median(count_frac, na.rm = TRUE),
    y_lo  = quantile(count_frac, 0.10, na.rm = TRUE),
    y_hi  = quantile(count_frac, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot(dist_sum, aes(x = depth, y = y_med)) +
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Depth distribution (median with 10-90% band; depth 0-100)",
    x = "Depth",
    y = "Fraction of bases"
  ) +
  theme_base

ggsave(file.path(out_dir, "depth_distribution_0_100.pdf"), p3, width = 6.8, height = 4.6)

chrom <- read_tsv(chrom_path, show_col_types = FALSE) %>%
  mutate(
    coverage = as.numeric(coverage),
    mean_depth = as.numeric(mean_depth),
    bin_mb = as.integer(bin_mb)
  ) %>%
  filter(is.finite(bin_mb)) %>%
  mutate(
    chrom = as.character(chrom),
    chrom = ifelse(str_detect(chrom, "^chr"), chrom, paste0("chr", chrom)),
    chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX")),
    cov_frac = coverage / 100
  ) %>%
  filter(!is.na(chrom), is.finite(cov_frac))

p4 <- ggplot(chrom, aes(x = bin_mb, y = cov_frac, group = sample_id)) +
  geom_line(linewidth = 0.25, alpha = 0.08) +
  facet_wrap(~ chrom, scales = "free_x", ncol = 6) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "Chromosome-wide 1 Mb coverage windows across samples",
    x = "1Mb genomic window",
    y = "Coverage"
  ) +
  theme_base +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "chrom_window_coverage_all_samples.pdf"), p4, width = 18, height = 10)
ggsave(file.path(out_dir, "chrom_window_coverage_all_samples.png"), p4, width = 18, height = 10, dpi = 300)

vaf <- read_tsv(vaf_path, show_col_types = FALSE) %>%
  mutate(vaf = as.numeric(vaf)) %>%
  filter(is.finite(vaf), vaf >= 0, vaf <= 1) %>%
  mutate(variant_class = factor(variant_class, levels = unique(variant_class)))

write_tsv(vaf, file.path(out_dir, "vaf_all_mutations_including_zero.tsv"))
vaf_exc0 <- vaf %>% filter(vaf > 0)
write_tsv(vaf_exc0, file.path(out_dir, "vaf_all_mutations_excluding_zero.tsv"))

plot_vaf_raw <- function(df, main_title, out_prefix) {
  p_txt <- NA_character_
  classes <- levels(droplevels(df$variant_class))
  if (length(classes) == 2) {
    tmp1 <- df %>% filter(variant_class == classes[1]) %>% pull(vaf)
    tmp2 <- df %>% filter(variant_class == classes[2]) %>% pull(vaf)
  } else {
    tmp1 <- numeric()
    tmp2 <- numeric()
  }
  if (length(tmp1) >= 2 && length(tmp2) >= 2) {
    pval <- suppressWarnings(wilcox.test(tmp1, tmp2, paired = FALSE)$p.value)
    p_txt <- paste0(classes[1], " versus ", classes[2], " Wilcoxon p = ", format(pval, digits = 3, scientific = TRUE))
  }

  p <- ggplot(df, aes(x = variant_class, y = vaf, fill = variant_class)) +
    geom_violin(trim = FALSE, width = 0.9, alpha = 0.35, color = "black", linewidth = 0.3) +
    ggbeeswarm::geom_quasirandom(
      width = 0.32,
      varwidth = TRUE,
      groupOnX = TRUE,
      size = 0.45,
      alpha = 0.22,
      stroke = 0
    ) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 7, color = "black") +
    labs(title = main_title, subtitle = p_txt, x = NULL, y = "Observed allele fraction") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_base +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), p, width = 6.5, height = 4.8)
  ggsave(file.path(out_dir, paste0(out_prefix, ".png")), p, width = 6.5, height = 4.8, dpi = 300)
}

plot_vaf_raw(vaf, "Observed allele-fraction distribution at supplied sites, including zero", "vaf_raw_including_zero")
plot_vaf_raw(vaf_exc0, "Observed allele-fraction distribution at supplied sites, excluding zero", "vaf_raw_excluding_zero")
