# prepare dataset
dataset <- generate_scenario(
  sim_num = 1L,
  IcMc_size = 10L, ItMc_size = 10L, IcMt_size = 10L, ItMt_size = 10L
)

multiple_dataset <- generate_scenario(
  sim_num = 3L,
  IcMc_size = 10L, ItMc_size = 10L, IcMt_size = 10L, ItMt_size = 10L
)

test_that("estimateEM_bootstrap_ci works", {
  boot_ci <- estimateEM_bootstrap_ci(dataset, boot_num = 5L, alpha = 0.05)
  expect_equal(nrow(boot_ci), 2L)
  expect_equal(ncol(boot_ci), 9L)
  expect_equal(boot_ci$percentile, c(0.025, 0.975))
})

test_that("estimate_ci_sequentially works", {
  boot_ci <- estimate_ci_sequentially(multiple_dataset, boot_num = 5L, alpha = 0.05, verbose = FALSE)
  expect_equal(nrow(boot_ci), 2L * 3L)
  expect_equal(ncol(boot_ci), 10L)
  expect_equal(boot_ci$percentile, rep(c(0.025, 0.975), 3L))
})

test_that("coverage works", {
  boot_ci <- estimate_ci_sequentially(multiple_dataset, boot_num = 5L, alpha = 0.05, verbose = FALSE)
  boot_coverage <- coverage(boot_ci)
  expect_equal(nrow(boot_coverage), 3L)
  expect_equal(ncol(boot_coverage), 13L)
})
