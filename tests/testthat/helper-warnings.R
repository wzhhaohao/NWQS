expect_small_boot_warning <- function(object) {
  pattern <- "'n_boot' is quite small; bootstrap percentile CI may be unstable\\."
  seen <- FALSE

  value <- withCallingHandlers(
    object,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        seen <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  expect_true(seen)
  value
}
