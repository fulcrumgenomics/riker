#!/usr/bin/env Rscript
# Curated performance figures embedded in the top-level README.
#
# These are presentation figures (Fulcrum-branded), distinct from the pipeline's
# diagnostic plot.R. They read a bench_summary.tsv produced by the benchmark
# pipeline and emit one PNG per figure.
#
# Usage (run from the repository root):
#   Rscript benchmark-pipeline/workflow/scripts/plot_readme_figures.R [bench_summary.tsv] [outdir]
#
# Defaults read the pipeline's merged results and write the PNGs the README embeds.

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
outdir <- if (length(args) >= 2) args[[2]] else "docs/images"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- Fulcrum palette (FG_Style_Guide.pdf) ------------------------------------
FG_GREEN <- "#38b44a" # riker — the hero series, everywhere
FG_NAVY  <- "#160052" # ink / Picard
FG_BLUE  <- "#26a8e0" # RustQC
FG_BLACK <- "#000000" # mosdepth (matches its logo's dominant color)
FG_GRAY  <- "#a5a8a7" # muted / gridlines
FONT     <- "Helvetica"

theme_fg <- function() {
  theme_minimal(base_family = FONT, base_size = 13) +
    theme(
      text             = element_text(color = FG_NAVY),
      plot.title       = element_text(face = "bold", size = 17, color = FG_NAVY),
      plot.subtitle    = element_text(size = 12, color = FG_NAVY, margin = margin(b = 10)),
      plot.caption     = element_text(size = 9, color = FG_GRAY, hjust = 0),
      axis.title       = element_text(size = 12, color = FG_NAVY),
      axis.text        = element_text(color = FG_NAVY),
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
speedup <- function(samp, prof, slow, fast = "riker") wall(samp, prof, slow) / wall(samp, prof, fast)

fmt_x      <- function(x) sprintf("%.1f×", x)
fmt_faster <- function(x) sprintf("%.1f× faster", x)

## ============================================================================
## Figure 1 — riker vs Picard, single-thread, tool-for-tool (speedup only)
## ============================================================================
wgs_cov <- c("HG02675_4x", "HG02675_15x", "HG02675_20x", "HG00188_30x", "HG02675_30x")
exome   <- c("HG002_av5", "HG003_av5", "HG004_av5")

multi_speed <- function(samples, prof, picard_a, picard_b) {
  mean(sapply(samples, function(x)
    (wall(x, prof, picard_a) + wall(x, prof, picard_b)) / wall(x, prof, "riker")))
}

lvl <- c("WGS", "Exome", "RNA · PE", "WGS · multi", "Exome · multi")
f1 <- tibble(
  label = factor(lvl, levels = lvl),
  speed = c(
    mean(sapply(wgs_cov, function(x) speedup(x, "wgs-t1", "picard-cwm"))),
    mean(sapply(exome,   function(x) speedup(x, "hybcap-only", "picard-chsm"))),
    speedup("ENCFF482AEE", "rna-t1", "picard-crsm"),
    multi_speed(wgs_cov, "wgs-bundle",    "picard-cmm", "picard-cwm"),
    multi_speed(exome,   "hybcap-bundle", "picard-cmm", "picard-chsm")
  )
)

p1 <- ggplot(f1, aes(label, speed)) +
  geom_col(fill = FG_GREEN, width = 0.68) +
  geom_text(aes(label = fmt_x(speed)), vjust = -0.45, family = FONT,
            fontface = "bold", size = 5, color = FG_NAVY) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), breaks = pretty_breaks(5)) +
  labs(
    title = "riker is 12–38× faster than Picard",
    subtitle = "Wall-clock speedup vs Picard — single tools at 1 thread, 'multi' runs at 4 threads",
    x = NULL, y = "× faster than Picard",
    caption = "single tools vs CollectWgsMetrics / CollectHsMetrics / CollectRnaSeqMetrics · multi vs CollectMultipleMetrics + the matching per-assay tool"
  ) +
  theme_fg() +
  theme(panel.grid.major.x = element_blank())

save_png(p1, "fig1_vs_picard.png", 8.8, 4.8)

## ============================================================================
## Figure 2a — WGS: riker vs mosdepth across coverage (single-thread)
## ============================================================================
cov_map <- c("4×" = "HG02675_4x", "15×" = "HG02675_15x",
             "20×" = "HG02675_20x", "30×" = "HG02675_30x")

f2a <- lapply(names(cov_map), function(cl) {
  samp <- cov_map[[cl]]
  tibble(
    cov = cl,
    tool = c("riker", "mosdepth"),
    secs = c(wall(samp, "wgs-t1", "riker"), wall(samp, "wgs-t1", "mosdepth"))
  )
}) %>% bind_rows() %>%
  mutate(cov = factor(cov, levels = names(cov_map)),
         tool = factor(tool, levels = c("mosdepth", "riker")))

f2a_lab <- f2a %>% group_by(cov) %>%
  summarise(top = max(secs),
            spd = secs[tool == "mosdepth"] / secs[tool == "riker"], .groups = "drop")

dodge <- position_dodge(width = 0.7)
p2a <- ggplot(f2a, aes(cov, secs, fill = tool)) +
  geom_col(position = dodge, width = 0.66) +
  geom_text(aes(label = sprintf("%.0f s", secs)), position = dodge, vjust = -0.4,
            family = FONT, size = 3.4, color = FG_NAVY) +
  geom_text(data = f2a_lab, aes(cov, top, label = fmt_faster(spd)), inherit.aes = FALSE,
            vjust = -2.1, family = FONT, fontface = "bold", size = 4.0, color = FG_GREEN) +
  scale_fill_manual(values = c(mosdepth = FG_BLACK, riker = FG_GREEN)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23))) +
  labs(title = "WGS: riker vs. mosdepth on 1000 Genomes HG02675",
       subtitle = "Single-threaded on 4×–30× BAMs",
       x = "Coverage", y = "Runtime (seconds)") +
  theme_fg() +
  theme(panel.grid.major.x = element_blank())

save_png(p2a, "fig2a_mosdepth_coverage.png", 7.2, 4.6)

## ============================================================================
## Figure 2b — WGS 30x: riker vs mosdepth across threads
## ============================================================================
thr_profs <- c("1" = "wgs-t1", "2" = "wgs-t2", "3" = "wgs-t3", "4" = "wgs-t4")
samp30 <- "HG02675_30x"

f2b <- lapply(names(thr_profs), function(tl) {
  prof <- thr_profs[[tl]]
  tibble(
    threads = tl,
    tool = c("riker", "mosdepth"),
    secs = c(wall(samp30, prof, "riker"), wall(samp30, prof, "mosdepth"))
  )
}) %>% bind_rows() %>%
  mutate(threads = factor(threads, levels = names(thr_profs)),
         tool = factor(tool, levels = c("mosdepth", "riker")))

f2b_lab <- f2b %>% group_by(threads) %>%
  summarise(top = max(secs),
            spd = secs[tool == "mosdepth"] / secs[tool == "riker"], .groups = "drop")

p2b <- ggplot(f2b, aes(threads, secs, fill = tool)) +
  geom_col(position = dodge, width = 0.66) +
  geom_text(aes(label = sprintf("%.0f s", secs)), position = dodge, vjust = -0.4,
            family = FONT, size = 3.4, color = FG_NAVY) +
  geom_text(data = f2b_lab, aes(threads, top, label = fmt_faster(spd)), inherit.aes = FALSE,
            vjust = -2.1, family = FONT, fontface = "bold", size = 4.0, color = FG_GREEN) +
  scale_fill_manual(values = c(mosdepth = FG_BLACK, riker = FG_GREEN)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.23))) +
  labs(title = "WGS: riker vs. mosdepth on 1000 Genomes HG02675",
       subtitle = "30× BAM across 1–4 threads",
       x = "Threads", y = "Runtime (seconds)") +
  theme_fg() +
  theme(panel.grid.major.x = element_blank())

save_png(p2b, "fig2b_mosdepth_threads.png", 7.2, 4.6)

## ============================================================================
## Figure 3 — RNA: riker vs Picard + RustQC (log runtime)
## ============================================================================
# The PE · 4 threads panel has no real Picard datum (CollectRnaSeqMetrics is
# single-threaded). We add a faded "ghost" Picard bar at its 1-thread runtime so
# readers see where it sits and that extra threads don't help it.
rna_rows <- tribble(
  ~grp,            ~grp_ord, ~tool,    ~samp,          ~prof,    ~ghost,
  "SE · 1 thread", 1, "Picard", "ENCFF028IUE", "rna-t1", FALSE,
  "SE · 1 thread", 1, "RustQC", "ENCFF028IUE", "rna-t1", FALSE,
  "SE · 1 thread", 1, "riker",  "ENCFF028IUE", "rna-t1", FALSE,
  "PE · 1 thread", 2, "Picard", "ENCFF482AEE", "rna-t1", FALSE,
  "PE · 1 thread", 2, "RustQC", "ENCFF482AEE", "rna-t1", FALSE,
  "PE · 1 thread", 2, "riker",  "ENCFF482AEE", "rna-t1", FALSE,
  "PE · 4 threads", 3, "Picard", "ENCFF482AEE", "rna-t1", TRUE,
  "PE · 4 threads", 3, "RustQC", "ENCFF482AEE", "rna-t4", FALSE,
  "PE · 4 threads", 3, "riker",  "ENCFF482AEE", "rna-t4", FALSE
)

f3 <- rna_rows %>%
  rowwise() %>%
  mutate(picard_tool = if (tool == "Picard") "picard-crsm" else if (tool == "RustQC") "rustqc" else "riker",
         secs = wall(samp, prof, picard_tool)) %>%
  ungroup() %>%
  group_by(grp) %>%
  mutate(spd = secs / secs[tool == "riker"]) %>%
  ungroup() %>%
  mutate(
    grp = fct_reorder(grp, grp_ord, .desc = TRUE), # .desc so SE ends up on top after coord_flip
    tool = factor(tool, levels = c("Picard", "RustQC", "riker")),
    lab = case_when(
      ghost           ~ sprintf("%.0f s · no multithreading", secs),
      tool == "riker" ~ sprintf("%.0f s · 1×", secs),
      TRUE            ~ sprintf("%.0f s · %.0f× slower", secs, spd)
    )
  )

dodge3 <- position_dodge2(preserve = "single", width = 0.75)
# group = tool pins both layers to the same dodge order; without it the alpha
# aesthetic on geom_col changes its grouping and the text labels misalign.
p3 <- ggplot(f3, aes(grp, secs, fill = tool, group = tool)) +
  geom_col(aes(alpha = ghost), position = dodge3, width = 0.72) +
  geom_text(aes(label = lab), position = dodge3,
            hjust = -0.06, family = FONT, size = 3.5, color = FG_NAVY) +
  scale_fill_manual(values = c(Picard = FG_NAVY, RustQC = FG_BLUE, riker = FG_GREEN)) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.25), guide = "none") +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.36)),
               breaks = c(0, 25, 100, 300, 600, 1000, 1400),
               minor_breaks = seq(100, 1400, 100),
               labels = function(x) paste0(x, " s")) +
  coord_flip() +
  labs(title = "RNA-seq: riker is 10–38× faster than Picard & RustQC",
       subtitle = "Wall-clock runtime on a √ scale (shorter is better)",
       x = NULL, y = "Runtime (seconds, √ scale)",
       caption = "Picard CollectRnaSeqMetrics is single-threaded; its faded 4-thread bar marks that same 1-thread runtime for reference.") +
  theme_fg() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(color = "grey85", size = 0.3))

save_png(p3, "fig3_rna.png", 8.6, 5.0)

message("done")
