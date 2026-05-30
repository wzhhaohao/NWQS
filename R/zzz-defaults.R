# Centralized default-value table for the v0.2.0 release. Every formal
# default on nwqs() and nwqs_boot() reads from NWQS_DEFAULTS rather than
# embedding a literal so the same value cannot drift between the two
# function signatures.
#
# To change a published default, update one entry here, run
# devtools::document(), and re-capture the relevant golden snapshots.

NWQS_DEFAULTS <- list(
  q                        = 4,
  df_spline                = 3,
  transform_type           = "percentile_rank",
  ties                     = "average",
  train_prop               = 0.6,
  rh                       = 10,
  n_permutation            = 30,
  n_boot                   = 100,
  rh_inner                 = 1,
  conf_level               = 0.95,
  seed                     = 1234,
  min_shape_sd             = 1e-8,
  zero_weight_action       = "na",
  contrast_target_pr       = 0.75,
  contrast_ref_pr          = 0.5,
  contrast_ref_q_bin       = 0,
  contrast_points_print_pr = c(0.25, 0.75, 0.95),
  effect_curve_grid        = seq(0, 1, by = 0.01),
  label_style              = "auto"
)
