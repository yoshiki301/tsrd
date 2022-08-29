# prepare dataset
X <- generate_scenario(sim_num = 1L)
X_IcMc <- subset(
  X, induction == "Ic" & maintenance == "Mc",
  select = c(time, cens, t_judge, is_induction)
)
X_IcMt <- subset(
  X, induction == "Ic" & maintenance == "Mt",
  select = c(time, cens, t_judge, is_induction)
)
X_censored <- generate_scenario(sim_num = 1L, censor_rate = 0.3)

# estimate results
pi_IcMc <- calc_responsibility(
  X_IM = X_IcMc,
  theta = 0.5, lambda_I_nr = 0.07,
  lambda_IM_r = 0.07, lambda_IM_nr = 0.10
)
pi_IcMt <- calc_responsibility(
  X_IM = X_IcMt,
  theta = 0.5, lambda_I_nr = 0.07,
  lambda_IM_r = 0.07, lambda_IM_nr = 0.10
)
theta_Ic <- update_theta(
  X_IMc = X_IcMc, X_IMt = X_IcMt,
  pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
)
lambda_Ic_nr <<- update_lambda_I(
  X_IMc = X_IcMc, X_IMt = X_IcMt,
  pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
)
lambda_IcMc <- update_lambda_IM(
  X_IM = X_IcMc, pi_IM = pi_IcMc
)

# calculate log-likelihood
loglikelihood_stage1 <- calc_loglikelihood_stage1(
  X_IMc = X_IcMc, X_IMt = X_IcMt,
  theta = theta_Ic, lambda_I_nr = lambda_Ic_nr
)
loglikelihood_stage2 <- calc_loglikelihood_stage2(
  X_IM = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
  lambda_IM_r = lambda_IcMc$r, lambda_IM_nr = lambda_IcMc$nr
)

test_that("calc_responsibility works", {
  expect_type(pi_IcMc, "double")
  expect_length(pi_IcMc, nrow(X_IcMc))
  expect_equal(all(0 <= pi_IcMc & pi_IcMc <= 1), T)
  expect_equal(sum(is.nan(pi_IcMc)) == 0, T)
})

test_that("update_theta works", {
  expect_type(theta_Ic, "double")
  expect_length(theta_Ic, 1L)
  expect_gte(theta_Ic, 0.0)
  expect_lte(theta_Ic, 1.0)
  expect_equal(!is.nan(theta_Ic), T)
})

test_that("update_lambda_I works", {
  expect_type(lambda_Ic_nr, "double")
  expect_length(lambda_Ic_nr, 1L)
  expect_gt(lambda_Ic_nr, 0.0)
  expect_equal(!is.nan(lambda_Ic_nr), T)
})

test_that("update_lambda_IM works", {
  expect_type(lambda_IcMc, "list")
  expect_type(lambda_IcMc$r, "double")
  expect_length(lambda_IcMc$r, 1L)
  expect_gt(lambda_IcMc$r, 0.0)
  expect_equal(!is.nan(lambda_IcMc$r) & !is.nan(lambda_IcMc$nr), T)
})

test_that("calc_loglikelihood_stage1 works", {
  expect_type(loglikelihood_stage1, "double")
  expect_length(loglikelihood_stage1, 1L)
  expect_equal(!is.nan(loglikelihood_stage1), T)
})

test_that("calc_loglikelihood_stage2 works", {
  expect_type(loglikelihood_stage2, "double")
  expect_length(loglikelihood_stage2, 1L)
  expect_equal(!is.nan(loglikelihood_stage2), T)
})

test_that("estimateEM works", {
  estimator <- estimateEM(dataset = X, max_iter = 100L)
  expect_type(estimator, "list")
  expect_equal(all(diff(estimator$loglikelihood) >= 0.0), T) # monotonous
})

test_that("estimateEM for censored data works", {
  estimator <- estimateEM(dataset = X_censored, max_iter = 100L)
  expect_type(estimator, "list")
  expect_equal(all(diff(estimator$loglikelihood) >= 0.0), T) # monotonous
})
