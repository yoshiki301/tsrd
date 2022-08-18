test_that("rexp_range works", {
  rand <- rexp_range(1000L, lambda = 1, range = c(0.5, 1))
  expect_equal(all(rand < 0.6931473), T)
})

test_that("sample_mixture_exp works", {
  rand <- sample_mixture_exp(1000L, theta = 0.5, lambda1 = 1, lambda2 = 2)
  expect_type(rand, "double")
  expect_length(rand, 1000L)
})

test_that("sample_tsrd_dataset works", {
  rand <- sample_tsrd_dataset(
    1000L, t_judge = 10, theta = 0.5,
    lambda_I_r = 1, lambda_I_nr = 1,
    lambda_IM_r = 2, lambda_IM_nr = 2
  )
  expect_type(rand, "double")
  expect_length(rand, 1000L)
})
