#!/usr/bin/env Rscript
# plot.R — produce 6 PDFs from bench.tsv + bench_summary.tsv.
#
# Usage:
#   Rscript plot.R <bench.tsv> <bench_summary.tsv> <plots_dir>

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript plot.R <bench.tsv> <bench_summary.tsv> <plots_dir>")
}
bench_path   <- args[1]
summary_path <- args[2]
plots_dir    <- args[3]
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

bench   <- read_tsv(bench_path,   show_col_types = FALSE)
summary <- read_tsv(summary_path, show_col_types = FALSE)

# Common style: small base size for dense panels, color palette mapping
# tool_family rather than tool token (so multiple picard variants share
# a hue but are distinguished by linetype/shape).
TOOL_FAMILIES <- c("riker", "picard", "mosdepth", "qualimap", "rustqc")
PALETTE <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"), TOOL_FAMILIES)

theme_bench <- function() {
  theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f4f4f4", color = NA))
}


# --- Plot 1: Wall-time leaderboard by profile ----------------------------
# Bars per (tool, sample), faceted by profile. Headline plot.

p1 <- summary %>%
  ggplot(aes(x = sample, y = wall_s_median, fill = tool_family)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           aes(group = tool)) +
  geom_errorbar(aes(ymin = wall_s_min, ymax = wall_s_max,
                    group = tool),
                position = position_dodge(0.8), width = 0.2,
                color = "grey30") +
  facet_wrap(~ profile, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = PALETTE) +
  labs(x = NULL, y = "wall time (s, median ± min/max across reps)",
       fill = "tool family",
       title = "Wall-time leaderboard by profile") +
  theme_bench() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(plots_dir, "01_wall_time_leaderboard.pdf"),
       p1, width = 10, height = 7)


# --- Plot 2: Throughput vs input size (carries the depth-series story) ---

p2 <- summary %>%
  ggplot(aes(x = sample_size_gb_median, y = throughput_mb_per_s_median,
             color = tool_family, shape = tool, group = tool)) +
  geom_point(size = 2.5) +
  geom_line(alpha = 0.6) +
  scale_x_log10(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  scale_color_manual(values = PALETTE) +
  facet_wrap(~ profile, scales = "free", ncol = 2) +
  labs(x = "input BAM size (GB, log scale)",
       y = "throughput (MB/s)",
       color = "tool family", shape = "tool",
       title = "Throughput vs input size",
       subtitle = "HG02675 4x/15x/20x/30x carries the size axis") +
  theme_bench()
ggsave(file.path(plots_dir, "02_throughput_vs_size.pdf"),
       p2, width = 10, height = 7)


# --- Plot 3: Max RSS vs input size (catches OOM-prone tools) -------------

p3 <- summary %>%
  ggplot(aes(x = sample_size_gb_median, y = max_rss_gb_median,
             color = tool_family, shape = tool, group = tool)) +
  geom_point(size = 2.5) +
  geom_line(alpha = 0.6) +
  scale_x_log10(labels = label_comma()) +
  scale_color_manual(values = PALETTE) +
  facet_wrap(~ profile, scales = "free", ncol = 2) +
  labs(x = "input BAM size (GB, log scale)",
       y = "peak RSS (GB)",
       color = "tool family", shape = "tool",
       title = "Peak RSS vs input size") +
  theme_bench()
ggsave(file.path(plots_dir, "03_max_rss_vs_size.pdf"),
       p3, width = 10, height = 7)


# --- Plot 4: Riker bundle vs sum-of-Picard ------------------------------
# For each (sample, bundle profile), compute:
#   wall_riker / sum(wall_picard-cmm + wall_picard-cwm/picard-chsm)
# Lower = riker faster than sum-of-Picard for the same coverage.

bundle_summary <- summary %>%
  filter(grepl("-bundle$", profile))

if (nrow(bundle_summary) > 0) {
  riker_bundle <- bundle_summary %>%
    filter(tool == "riker") %>%
    select(sample, profile, riker_wall = wall_s_median)

  picard_bundle <- bundle_summary %>%
    filter(tool_family == "picard") %>%
    group_by(sample, profile) %>%
    summarise(picard_total_wall = sum(wall_s_median, na.rm = TRUE),
              .groups = "drop")

  bundle_ratio <- riker_bundle %>%
    inner_join(picard_bundle, by = c("sample", "profile")) %>%
    mutate(ratio = riker_wall / picard_total_wall)

  p4 <- bundle_ratio %>%
    ggplot(aes(x = sample, y = ratio, fill = profile)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = sprintf("%.2fx", ratio)),
              position = position_dodge(0.8), vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("wgs-bundle" = "#0072B2",
                                  "hybcap-bundle" = "#E69F00")) +
    labs(x = NULL, y = "wall(riker bundle) / sum(wall(picard parts))",
         fill = "profile",
         title = "Riker single-pass vs sum of Picard tools",
         subtitle = "ratio < 1 = riker is faster than running the equivalent Picard tools sequentially") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(plots_dir, "04_riker_bundle_vs_picard.pdf"),
         p4, width = 10, height = 6)
} else {
  # Empty placeholder so snakemake's rule output exists.
  pdf(file.path(plots_dir, "04_riker_bundle_vs_picard.pdf")); plot.new(); dev.off()
}


# --- Plot 5: Riker multi vs Picard CMM head-to-head ---------------------
# Single-pass vs single-pass on the overlapping subset (the bundle
# profiles).

cmm_h2h <- summary %>%
  filter(grepl("-bundle$", profile),
         tool %in% c("riker", "picard-cmm")) %>%
  select(sample, profile, tool, wall_s_median, max_rss_gb_median, throughput_mb_per_s_median)

if (nrow(cmm_h2h) > 0) {
  p5 <- cmm_h2h %>%
    ggplot(aes(x = sample, y = wall_s_median, fill = tool)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    facet_wrap(~ profile, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("riker" = "#0072B2", "picard-cmm" = "#E69F00")) +
    labs(x = NULL, y = "wall time (s)",
         fill = "tool",
         title = "Single-pass head-to-head: riker multi vs Picard CollectMultipleMetrics",
         subtitle = "Picard CMM is run with USE_ASYNC_IO_*=true") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(plots_dir, "05_riker_vs_picard_cmm.pdf"),
         p5, width = 10, height = 6)
} else {
  pdf(file.path(plots_dir, "05_riker_vs_picard_cmm.pdf")); plot.new(); dev.off()
}


# --- Plot 6: Per-comparator ratio (riker = 1.0 baseline) ----------------
# wall_other / wall_riker per (sample, profile, tool). Compresses the
# whole matrix into one plot.

riker_ref <- summary %>%
  filter(tool == "riker") %>%
  select(sample, profile, riker_wall = wall_s_median)

ratio_df <- summary %>%
  filter(tool != "riker") %>%
  inner_join(riker_ref, by = c("sample", "profile")) %>%
  mutate(ratio = wall_s_median / riker_wall)

if (nrow(ratio_df) > 0) {
  p6 <- ratio_df %>%
    ggplot(aes(x = tool, y = ratio, fill = tool_family)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    facet_grid(profile ~ sample, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = PALETTE) +
    labs(x = NULL, y = "wall(comparator) / wall(riker)",
         fill = "tool family",
         title = "Per-comparator ratio (riker = 1.0)") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          strip.text.x = element_text(size = 7))
  ggsave(file.path(plots_dir, "06_per_comparator_ratio.pdf"),
         p6, width = 14, height = 8)
} else {
  pdf(file.path(plots_dir, "06_per_comparator_ratio.pdf")); plot.new(); dev.off()
}

cat("Wrote 6 PDFs to", plots_dir, "\n")
