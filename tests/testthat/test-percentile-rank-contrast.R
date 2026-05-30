test_that(".contrast_point_label formats percentile_rank points as P{round*100}", {
  expect_equal(NWQS:::.contrast_point_label(0,   "percentile_rank", "P"), "P0")
  expect_equal(NWQS:::.contrast_point_label(0.5, "percentile_rank", "P"), "P50")
  expect_equal(NWQS:::.contrast_point_label(1/3, "percentile_rank", "P"), "P33")
  expect_equal(NWQS:::.contrast_point_label(2/3, "percentile_rank", "P"), "P67")
  expect_equal(NWQS:::.contrast_point_label(1,   "percentile_rank", "P"), "P100")
})

test_that(".contrast_point_label preserves Q labels for q_bin", {
  expect_equal(NWQS:::.contrast_point_label(0, "q_bin", "Q"), "Q1")
  expect_equal(NWQS:::.contrast_point_label(3, "q_bin", "Q"), "Q4")
})

test_that(".contrast_point_label numeric style returns trimmed numeric", {
  expect_equal(NWQS:::.contrast_point_label(0.25, "percentile_rank", "numeric"), "0.25")
  expect_equal(NWQS:::.contrast_point_label(3,    "q_bin",            "numeric"), "3")
})

test_that(".contrast_pair_label joins with _vs_", {
  expect_equal(
    NWQS:::.contrast_pair_label(0.75, 0.25, "percentile_rank", "P"),
    "P75_vs_P25"
  )
  expect_equal(
    NWQS:::.contrast_pair_label(3, 0, "q_bin", "Q"),
    "Q4_vs_Q1"
  )
})

test_that(".label_style_default picks P for percentile_rank, Q for q_bin", {
  expect_equal(NWQS:::.label_style_default("percentile_rank"), "P")
  expect_equal(NWQS:::.label_style_default("q_bin"),           "Q")
})

test_that(".validate_pr_points rejects out-of-[0,1] in percentile_rank", {
  expect_error(NWQS:::.validate_pr_points(c(0.5, 1.2), "percentile_rank"), "\\[0, 1\\]")
  expect_silent(NWQS:::.validate_pr_points(c(0, 0.5, 1), "percentile_rank"))
  expect_silent(NWQS:::.validate_pr_points(c(0, 3),     "q_bin"))
})
