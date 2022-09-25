# prepare dataset
dataset <- generate_scenario(
  sim_num = 1L,
  IcMc_size = 10L, ItMc_size = 10L, IcMt_size = 10L, ItMt_size = 10L
)

test_that("calc_responsibility works", {
  boot_ci <- estimateEM_bootstrap_ci(dataset, boot_num = 5L, alpha = 0.05)
  expect_equal(nrow(boot_ci), 2)
  expect_equal(ncol(boot_ci), 5)
  expect_equal(boot_ci$percentile, c(0.025, 0.975))
})
