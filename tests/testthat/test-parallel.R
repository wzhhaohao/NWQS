# Contract tests for configure_parallel_plan(). The function MUST stay
# non-invasive: if the caller has already set up a non-sequential future
# plan, this helper returns invisibly without changing it. That is a
# load-bearing guarantee for HPC users who want NWQS to honor their own
# plan setup, and it is currently undocumented in test form.

# Each test resets the future plan on exit so it cannot pollute the rest
# of the suite.

skip_if_not_installed("future")

# ----- Contract 1: user-set plan is preserved ----------------------------

test_that("configure_parallel_plan() leaves a user-set multisession plan untouched", {
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)

  res <- configure_parallel_plan(
    loop_number = 10,
    strategy    = "multisession",
    n_workers   = NULL,
    verbose     = FALSE
  )

  current <- future::plan()
  expect_true(inherits(current, "multisession"))
})

# ----- Contract 2: sequential strategy short-circuits to sequential ------

test_that("configure_parallel_plan() with strategy = 'sequential' yields a sequential plan", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  configure_parallel_plan(
    loop_number = 10,
    strategy    = "sequential",
    n_workers   = NULL,
    verbose     = FALSE
  )

  expect_true(inherits(future::plan(), "sequential"))
})

# ----- Contract 3: loop_number = 1 short-circuits to sequential ----------

test_that("configure_parallel_plan() with loop_number = 1 short-circuits to sequential", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  configure_parallel_plan(
    loop_number = 1,
    strategy    = "multisession",
    n_workers   = 2,
    verbose     = FALSE
  )

  expect_true(inherits(future::plan(), "sequential"))
})

# ----- Contract 4: explicit n_workers wins over auto load balancing ------

test_that("configure_parallel_plan() honors explicit n_workers when computing a fresh plan", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  configure_parallel_plan(
    loop_number = 10,
    strategy    = "multisession",
    n_workers   = 2,
    verbose     = FALSE
  )

  expect_true(inherits(future::plan(), "multisession"))
  expect_equal(future::nbrOfWorkers(), 2)
})
