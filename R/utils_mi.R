#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Shared Functions: Rubin's Rules for Multiply Imputed DXA Data
#------------------------------------------------------------------------------------------

# NHANES releases 1999-2004 DXA body composition as five multiply imputed
# copies per participant. A model is fit once per copy and the results are
# combined with Rubin's rules: the pooled estimate is the mean, and the
# pooled variance adds the between-copy spread to the average within-copy
# variance.
pool_rubin <- function(estimates, variances) {
  m <- length(estimates)
  stopifnot(m == length(variances), m >= 2)
  q_bar <- mean(estimates)
  within <- mean(variances)
  between <- stats::var(estimates)
  total <- within + (1 + 1 / m) * between
  # Barnard-Rubin style degrees of freedom (large-sample version).
  r <- (1 + 1 / m) * between / within
  df <- if (between == 0) Inf else (m - 1) * (1 + 1 / r)^2
  t_crit <- stats::qt(0.975, df)
  data.frame(estimate = q_bar, se = sqrt(total),
             lower = q_bar - t_crit * sqrt(total), upper = q_bar + t_crit * sqrt(total),
             df = df)
}
