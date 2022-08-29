# prepare dataset
X <- generate_scenario(
  sim_num = 10L,
  IcMc_size = 20L,
  IcMt_size = 30L,
  ItMc_size = 40L,
  ItMt_size = 50L
)
num_cores <- 2L
if (parallel::detectCores() == 1) { num_cores <- 1 }

test_that("estimate_parallely works", {
  estimators <- estimate_parallely(
    dataset = X,
    num_cores = num_cores,
    verbose = T,
    max_iter = 100L
  )
  expect_type(estimators, "list")
  expect_equal(nrow(estimators), 10L)
})
