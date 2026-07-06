#!/usr/bin/env Rscript
# Denser performance figures for benchmark-pipeline/README.md. Same Fulcrum
# styling as the top-level README's plot_readme_figures.R, but these carry the
# full thread x coverage sweeps and the memory story rather than the headline
# slices.
#
# Usage (run from the repository root):
#   Rscript benchmark-pipeline/workflow/scripts/plot_benchmark_figures.R [bench_summary.tsv] [outdir]

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
summary_path <- if (length(args) >= 1) args[[1]] else "benchmark-pipeline/results/bench_summary.tsv"
outdir <- if (length(args) >= 2) args[[2]] else "benchmark-pipeline/docs/images"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- Fulcrum palette ---------------------------------------------------------
FG_GREEN <- "#38b44a" # riker
FG_NAVY  <- "#160052" # Picard / ink
FG_BLUE  <- "#26a8e0" # RustQC
FG_BLACK <- "#000000" # mosdepth
FONT     <- "Helvetica"

theme_fg <- function() {
  theme_minimal(base_family = FONT, base_size = 13) +
    theme(
      text             = element_text(color = FG_NAVY),
      plot.title       = element_text(face = "bold", size = 17, color = FG_NAVY),
      plot.subtitle    = element_text(size = 12, color = FG_NAVY, margin = margin(b = 10)),
      plot.caption     = element_text(size = 9, color = "#7a7a90", hjust = 0),
      axis.title       = element_text(size = 12, color = FG_NAVY),
      axis.text        = element_text(color = FG_NAVY),
      strip.text       = element_text(face = "bold", size = 12, color = FG_NAVY),
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 12),
      plot.margin      = margin(14, 18, 10, 12),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

save_png <- function(p, file, w, h) {
  path <- file.path(outdir, file)
  ragg::agg_png(path, width = w, height = h, units = "in", res = 200, background = "white")
  print(p)
  invisible(dev.off())
  message("wrote ", path)
}

s <- read_tsv(summary_path, show_col_types = FALSE)
wall <- function(samp, prof, tl) {
  v <- s %>% filter(sample == samp, profile == prof, tool == tl) %>% pull(wall_s_median)
  if (length(v) == 0) NA_real_ else as.numeric(v[[1]])
}

## ============================================================================
## Figure A — WGS thread scaling vs mosdepth, faceted by coverage
## ============================================================================
cov_map <- c("4× coverage" = "HG02675_4x", "15× coverage" = "HG02675_15x",
             "20× coverage" = "HG02675_20x", "30× coverage" = "HG02675_30x")
thr_map <- c("1" = "wgs-t1", "2" = "wgs-t2", "3" = "wgs-t3", "4" = "wgs-t4")

wgs <- expand_grid(cov = names(cov_map), thr = names(thr_map),
                   tool = c("riker", "mosdepth")) %>%
  rowwise() %>%
  mutate(secs = wall(cov_map[[cov]], thr_map[[thr]], tool)) %>%
  ungroup() %>%
  mutate(cov = factor(cov, levels = names(cov_map)),
         threads = as.integer(thr),
         tool = factor(tool, levels = c("mosdepth", "riker")))

pA <- ggplot(wgs, aes(threads, secs, color = tool)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.4) +
  facet_wrap(~cov, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c(mosdepth = FG_BLACK, riker = FG_GREEN)) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "WGS depth: riker vs. mosdepth, full thread × coverage sweep",
       subtitle = "1000 Genomes HG02675, cold cache, harmonized reads · riker is 1.2–1.5× faster at every point",
       x = "Threads (riker --threads N; mosdepth -t N-1)", y = "Runtime (seconds)") +
  theme_fg()

save_png(pA, "fig_wgs_scaling.png", 8.6, 6.2)

## ============================================================================
## Figure B — RNA thread scaling vs RustQC, Picard as single-threaded reference
## ============================================================================
lay_map <- c("single-end" = "ENCFF028IUE", "paired-end" = "ENCFF482AEE")
rthr    <- c("1" = "rna-t1", "2" = "rna-t2", "4" = "rna-t4")

rna <- expand_grid(layout = names(lay_map), thr = names(rthr),
                   tool = c("riker", "rustqc")) %>%
  rowwise() %>%
  mutate(secs = wall(lay_map[[layout]], rthr[[thr]], tool)) %>%
  ungroup() %>%
  mutate(layout = factor(layout, levels = names(lay_map)),
         threads = as.integer(thr),
         tool = factor(recode(tool, rustqc = "RustQC"), levels = c("riker", "RustQC")))

picard_ref <- tibble(layout = factor(names(lay_map), levels = names(lay_map)),
                     secs = sapply(lay_map, function(x) wall(x, "rna-t1", "picard-crsm")))

pB <- ggplot(rna, aes(threads, secs, color = tool)) +
  geom_hline(data = picard_ref, aes(yintercept = secs),
             linetype = "dashed", color = FG_NAVY, size = 0.5) +
  geom_text(data = picard_ref, aes(x = 4, y = secs, label = "Picard (single-threaded)"),
            inherit.aes = FALSE, vjust = -0.6, hjust = 1, family = FONT,
            size = 3.3, color = FG_NAVY) +
  geom_line(size = 1.1) +
  geom_point(size = 2.4) +
  facet_wrap(~layout, nrow = 1) +
  scale_color_manual(values = c(riker = FG_GREEN, RustQC = FG_BLUE)) +
  scale_x_continuous(breaks = c(1, 2, 4)) +
  scale_y_log10(labels = function(x) paste0(x, " s"),
                breaks = c(10, 30, 100, 300, 1000)) +
  labs(title = "RNA-seq: riker vs. RustQC thread scaling",
       subtitle = "ENCODE polyA BAMs · log runtime · Picard CollectRnaSeqMetrics cannot thread (dashed reference)",
       x = "Threads", y = "Runtime (log scale)") +
  theme_fg()

save_png(pB, "fig_rna_scaling.png", 8.6, 4.8)

## ============================================================================
## Figure C — peak memory footprint by tool (across the whole benchmark)
## ============================================================================
fam_label <- c(riker = "riker", mosdepth = "mosdepth", picard = "Picard", rustqc = "RustQC")
mem <- s %>%
  filter(tool_family %in% names(fam_label)) %>%
  group_by(tool_family) %>%
  summarise(peak = max(as.numeric(max_rss_gb_median), na.rm = TRUE), .groups = "drop") %>%
  mutate(tool = factor(fam_label[tool_family],
                       levels = c("riker", "mosdepth", "Picard", "RustQC")))

pC <- ggplot(mem, aes(tool, peak, fill = tool)) +
  geom_col(width = 0.64) +
  geom_text(aes(label = sprintf("%.1f GB", peak)), vjust = -0.4,
            family = FONT, fontface = "bold", size = 4.6, color = FG_NAVY) +
  scale_fill_manual(values = c(riker = FG_GREEN, mosdepth = FG_BLACK,
                               Picard = FG_NAVY, RustQC = FG_BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Peak memory footprint",
       subtitle = "Largest resident set size reached by each tool anywhere in the benchmark (lower is better)",
       x = NULL, y = "Peak RSS (GB)") +
  theme_fg() +
  theme(panel.grid.major.x = element_blank())

save_png(pC, "fig_memory.png", 8.0, 4.6)

message("done")
