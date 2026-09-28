#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Figure 1 (Gap by Measure) and Figure 2 (Per-Marker Decomposition)
#------------------------------------------------------------------------------------------

# Figure 1: Black-White gap in each aging measure, NHANES vs. Graf et al.
# Figure 2: the blood PhenoAge gap split by marker.
# Colors are the reference categorical slots 1-2 and the blue/red diverging
# pair, checked with a colorblind-separation validator (all checks pass).

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

SURFACE <- "#fcfcfb"
INK <- "#0b0b0b"
INK_2 <- "#52514e"
GRID <- "#e6e5e1"
BLUE <- "#2a78d6"
ORANGE <- "#eb6834"
RED <- "#e34948"

theme_repo <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = SURFACE, colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
      axis.text = element_text(colour = INK_2),
      axis.title = element_text(colour = INK_2),
      plot.title = element_text(colour = INK, face = "bold"),
      plot.subtitle = element_text(colour = INK_2),
      plot.caption = element_text(colour = INK_2, hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      legend.position = "top", legend.justification = "left",
      legend.text = element_text(colour = INK_2), legend.title = element_blank()
    )
}
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Figure 1 -----------------------------------------------------------------
g <- read.csv("output/tables/gap_by_clock.csv")
order <- c("PhenoAge (blood chemistry)", "PhenoAge (epigenetic)", "GrimAge", "GrimAge2",
           "DunedinPoAm", "Horvath clock", "Hannum clock")
fig1_data <- bind_rows(
  g |> transmute(label, d, lo = d_lower, hi = d_upper, source = "NHANES 1999-2002 (this repo)"),
  g |> filter(!is.na(graf_cohens_d)) |>
    transmute(label, d = graf_cohens_d, lo = graf_d_lower, hi = graf_d_upper,
              source = "HRS 2016 (Graf et al. 2022)")
) |>
  mutate(label = factor(label, levels = rev(order)),
         source = factor(source, levels = c("NHANES 1999-2002 (this repo)", "HRS 2016 (Graf et al. 2022)")))

fig1 <- ggplot(fig1_data, aes(d, label, colour = source, shape = source)) +
  geom_vline(xintercept = 0, colour = INK_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(BLUE, ORANGE)) +
  scale_shape_manual(values = c(16, 17)) +
  annotate("text", x = 0.02, y = 7.45, label = "Black participants measure older →",
           hjust = 0, size = 3, colour = INK_2) +
  annotate("text", x = -0.02, y = 7.45, label = "← younger",
           hjust = 1, size = 3, colour = INK_2) +
  labs(
    title = "Same people, two PhenoAges, two answers",
    subtitle = "Black-White difference in each measure, adjusted for age (Cohen's d, 95% CI)",
    x = "Standardized difference (Cohen's d)", y = NULL,
    caption = "NHANES: n = 1,556 adults 50+ with both blood chemistry and DNA methylation, survey-weighted.\nHRS values from Graf et al. 2022, Web Table 9. GrimAge2 was not in Graf's table."
  ) +
  theme_repo()
ggsave("output/figures/fig1_gap_by_clock.png", fig1, width = 8, height = 5, dpi = 200, bg = SURFACE)

# ---- Figure 2 -----------------------------------------------------------------
dec <- read.csv("output/tables/blood_phenoage_decomposition.csv") |>
  mutate(label = factor(label, levels = rev(label)),
         direction = if_else(contribution_years > 0, "Makes Black participants look older",
                             "Makes Black participants look younger"),
         text = sprintf("%+.2f", contribution_years),
         hjust = if_else(contribution_years > 0, -0.15, 1.15))
total <- sum(dec$contribution_years)

fig2 <- ggplot(dec, aes(contribution_years, label, fill = direction)) +
  geom_vline(xintercept = 0, colour = INK_2, linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_text(aes(label = text, hjust = hjust), size = 3.2, colour = INK) +
  scale_fill_manual(values = c(RED, BLUE)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(
    title = sprintf("Where the %.1f-year blood PhenoAge gap comes from", total),
    subtitle = "Each marker's contribution, in years; contributions add up exactly to the total",
    x = "Contribution to the Black-White gap (years)", y = NULL,
    caption = "Contribution = marker's age-adjusted Black-White difference x its weight in Levine's formula (converted to years)."
  ) +
  theme_repo()
ggsave("output/figures/fig2_decomposition.png", fig2, width = 8, height = 5, dpi = 200, bg = SURFACE)
