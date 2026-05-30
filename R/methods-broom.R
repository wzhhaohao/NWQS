# broom::tidy and broom::glance support for nwqs and nwqs_boot objects.
#
# Registered conditionally in R/zzz.R when the user has the `broom`
# (or `generics`) package installed, so NWQS does not need a hard
# dependency on broom for users who do not want it.
#
# Outputs are plain data.frames; users on the tidyverse stack can wrap
# with tibble::as_tibble() if they want a true tibble.

tidy_nwqs <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  fit_coefs <- x$fit$coefficients

  if (is.data.frame(fit_coefs)) {
    coef_df <- fit_coefs
  } else {
    coef_df <- as.data.frame(fit_coefs)
  }

  z_col <- if ("z value" %in% names(coef_df)) "z value"
           else if ("t value" %in% names(coef_df)) "t value"
           else NA_character_
  p_col <- grep("^Pr\\(", names(coef_df), value = TRUE)
  if (length(p_col) == 0) p_col <- NA_character_ else p_col <- p_col[1]

  out <- data.frame(
    term      = rownames(coef_df),
    estimate  = coef_df$Estimate,
    std.error = coef_df[["Std. Error"]],
    statistic = if (!is.na(z_col)) coef_df[[z_col]] else NA_real_,
    p.value   = if (!is.na(p_col)) coef_df[[p_col]] else NA_real_,
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    ci <- tryCatch(
      suppressWarnings(confint.nwqs(x, level = conf.level)),
      error = function(e) NULL
    )
    if (!is.null(ci)) {
      out$conf.low  <- ci[out$term, 1]
      out$conf.high <- ci[out$term, 2]
    }
  }

  rownames(out) <- NULL
  out
}

glance_nwqs <- function(x, ...) {
  fit_metrics <- x$fit
  data.frame(
    n              = nrow(x$data),
    family         = x$family,
    transform_type = if (!is.null(x$transform_type)) x$transform_type else NA_character_,
    rh             = x$rh,
    aic            = if (!is.null(fit_metrics$aic)) fit_metrics$aic else NA_real_,
    deviance       = if (!is.null(fit_metrics$deviance)) fit_metrics$deviance else NA_real_,
    null.deviance  = if (!is.null(fit_metrics$null.deviance)) fit_metrics$null.deviance else NA_real_,
    df.residual    = if (!is.null(fit_metrics$df.residual)) fit_metrics$df.residual else NA_real_,
    stringsAsFactors = FALSE
  )
}

tidy_nwqs_boot <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  if (is.null(conf.level)) {
    conf.level <- if (!is.null(x$conf_level)) x$conf_level else 0.95
  }
  coef_names <- names(x$mean_coefs)

  vmat <- tryCatch(
    vcov.nwqs_boot(x),
    error = function(e) NULL
  )
  std_errors <- if (!is.null(vmat)) sqrt(diag(vmat))[coef_names] else NA_real_

  out <- data.frame(
    term      = coef_names,
    estimate  = unname(x$mean_coefs),
    std.error = unname(std_errors),
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    ci <- tryCatch(
      confint.nwqs_boot(x, level = conf.level),
      error = function(e) NULL
    )
    if (!is.null(ci)) {
      out$conf.low  <- ci[out$term, 1]
      out$conf.high <- ci[out$term, 2]
    }
  }

  rownames(out) <- NULL
  out
}

glance_nwqs_boot <- function(x, ...) {
  data.frame(
    n              = nrow(x$data),
    family         = x$family,
    transform_type = if (!is.null(x$transform_type)) x$transform_type else NA_character_,
    n_boot         = x$n_boot,
    n_success      = x$n_success,
    conf_level     = if (!is.null(x$conf_level)) x$conf_level else NA_real_,
    rh_inner       = if (!is.null(x$rh_inner)) x$rh_inner else NA_integer_,
    stringsAsFactors = FALSE
  )
}
