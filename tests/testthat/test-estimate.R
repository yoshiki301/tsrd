# prepare dataset
X <- generate_scenario(
  sim_num = 1L,
  IcMc_size = 100L,
  ItMc_size = 200L,
  IcMt_size = 300L,
  ItMt_size = 400L
)
X_IcMc <- subset(
  X, induction == "Ic" & maintenance == "Mc",
  select = c(time, cens, t_judge, is_induction)
)
X_IcMt <- subset(
  X, induction == "Ic" & maintenance == "Mt",
  select = c(time, cens, t_judge, is_induction)
)
X_ItMc <- subset(
  X, induction == "It" & maintenance == "Mc",
  select = c(time, cens, t_judge, is_induction)
)
X_ItMt <- subset(
  X, induction == "It" & maintenance == "Mt",
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
pi_ItMc <- calc_responsibility(
  X_IM = X_ItMc,
  theta = 0.5, lambda_I_nr = 0.07,
  lambda_IM_r = 0.07, lambda_IM_nr = 0.10
)
pi_ItMt <- calc_responsibility(
  X_IM = X_ItMt,
  theta = 0.5, lambda_I_nr = 0.07,
  lambda_IM_r = 0.07, lambda_IM_nr = 0.10
)
theta_Ic <- update_theta(
  X_IMc = X_IcMc, X_IMt = X_IcMt,
  pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
)
theta_It <- update_theta(
  X_IMc = X_ItMc, X_IMt = X_ItMt,
  pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
)
lambda_Ic_nr <- update_lambda_I(
  X_IMc = X_IcMc, X_IMt = X_IcMt,
  pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
)
lambda_It_nr <- update_lambda_I(
  X_IMc = X_ItMc, X_IMt = X_ItMt,
  pi_IMc = pi_ItMc, pi_IMt = pi_ItMt
)
lambda_IcMc <- update_lambda_IM(
  X_IM = X_IcMc, pi_IM = pi_IcMc
)
lambda_IcMt <- update_lambda_IM(
  X_IM = X_IcMt, pi_IM = pi_IcMt
)
lambda_ItMc <- update_lambda_IM(
  X_IM = X_ItMc, pi_IM = pi_ItMc
)
lambda_ItMt <- update_lambda_IM(
  X_IM = X_ItMt, pi_IM = pi_ItMt
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

test_that("calc_Fisher_information works", {
  parameters <- list(
    pi_IcMc = pi_IcMc,
    pi_IcMt = pi_IcMt,
    pi_ItMc = pi_ItMc,
    pi_ItMt = pi_ItMt,
    theta_Ic = theta_Ic,
    theta_It = theta_It,
    lambda_Ic_r = 0.0,
    lambda_Ic_nr = lambda_Ic_nr,
    lambda_It_r = 0.0,
    lambda_It_nr = lambda_It_nr,
    lambda_IcMc_r = lambda_IcMc$r,
    lambda_IcMc_nr = lambda_IcMc$nr,
    lambda_ItMc_r = lambda_ItMc$r,
    lambda_ItMc_nr = lambda_ItMc$nr,
    lambda_IcMt_r = lambda_IcMt$r,
    lambda_IcMt_nr = lambda_IcMt$nr,
    lambda_ItMt_r = lambda_ItMt$r,
    lambda_ItMt_nr = lambda_ItMt$nr
  )
  fisher_information <- calc_Fisher_information(
    dataset = X, parameters = parameters
  )
  expect_type(fisher_information, "list")
  expect_equal(sum(is.na(fisher_information)), 0L)
  expect_equal(all(fisher_information > 0), T) # diagonal elements must have positivity
  # TODO: check non-diagonal elements
})

test_that("estimateEM works", {
  estimator <- estimateEM(dataset = X, max_iter = 100L)
  expect_type(estimator, "list")
  expect_equal(all(diff(estimator$loglikelihood) >= 0.0), T) # monotonous
  expected_names <- c(
    "step", "loglikelihood",
    "pi_IcMc", "pi_IcMt", "pi_ItMc", "pi_ItMt",
    "theta_Ic", "theta_It",
    "lambda_Ic_r", "lambda_Ic_nr",
    "lambda_It_r", "lambda_It_nr",
    "lambda_IcMc_r", "lambda_IcMc_nr",
    "lambda_ItMc_r", "lambda_ItMc_nr",
    "lambda_IcMt_r", "lambda_IcMt_nr",
    "lambda_ItMt_r", "lambda_ItMt_nr",
    "theta_Ic_info", "theta_It_info",
    "lambda_Ic_nr_info",
    "lambda_It_nr_info",
    "lambda_IcMc_r_info", "lambda_IcMc_nr_info",
    "lambda_ItMc_r_info", "lambda_ItMc_nr_info",
    "lambda_IcMt_r_info", "lambda_IcMt_nr_info",
    "lambda_ItMt_r_info", "lambda_ItMt_nr_info"
  )
  expect_equal(names(estimator), expected_names)
})

test_that("estimateEM for censored data works", {
  estimator <- estimateEM(dataset = X_censored, max_iter = 100L)
  expect_type(estimator, "list")
  expect_equal(all(diff(estimator$loglikelihood) >= 0.0), T) # monotonous
})

test_that("estimateEM to fix theta", {
  estimator <- estimateEM(dataset = X, max_iter = 100L, option = "fix_theta")
  expect_type(estimator, "list")
  expect_equal(all(diff(estimator$loglikelihood) >= 0.0), T) # monotonous
  expect_equal(estimator$theta_Ic[1], 0.65)
  expect_equal(estimator$theta_It[1], 0.75)
})

test_that("error in estimateEM", {
  expect_error(
    estimateEM(dataset = X, max_iter = 100L, option = "invalid_option_param"),
    "option parameter is invalid."
  )
})

test_that("estimateEM_as_frame works", {
  estimator <- estimateEM_as_frame(dataset = X, max_iter = 100L)
  expect_type(estimator, "list")
  expected_colnames <- c(
    "theta_Ic", "theta_It",
    "lambda_Ic_r", "lambda_Ic_nr",
    "lambda_It_r", "lambda_It_nr",
    "lambda_IcMc_r", "lambda_IcMc_nr",
    "lambda_ItMc_r", "lambda_ItMc_nr",
    "lambda_IcMt_r", "lambda_IcMt_nr",
    "lambda_ItMt_r", "lambda_ItMt_nr",
    "theta_Ic_info", "theta_It_info",
    "lambda_Ic_nr_info",
    "lambda_It_nr_info",
    "lambda_IcMc_r_info", "lambda_IcMc_nr_info",
    "lambda_ItMc_r_info", "lambda_ItMc_nr_info",
    "lambda_IcMt_r_info", "lambda_IcMt_nr_info",
    "lambda_ItMt_r_info", "lambda_ItMt_nr_info",
    "last_loglikelihood",
    "n_IcMc", "n_ItMc", "n_IcMt", "n_ItMt",
    "t_judge"
  )
  expect_equal(colnames(estimator), expected_colnames)
})

test_that("estimate_sequentially works", {
  estimators <- estimate_sequentially(dataset = X, verbose = FALSE, max_iter = 100L)
  expect_type(estimators, "list")
  expected_colnames <- c(
    "id",
    "theta_Ic", "theta_It",
    "lambda_Ic_r", "lambda_Ic_nr",
    "lambda_It_r", "lambda_It_nr",
    "lambda_IcMc_r", "lambda_IcMc_nr",
    "lambda_ItMc_r", "lambda_ItMc_nr",
    "lambda_IcMt_r", "lambda_IcMt_nr",
    "lambda_ItMt_r", "lambda_ItMt_nr",
    "theta_Ic_info", "theta_It_info",
    "lambda_Ic_nr_info",
    "lambda_It_nr_info",
    "lambda_IcMc_r_info", "lambda_IcMc_nr_info",
    "lambda_ItMc_r_info", "lambda_ItMc_nr_info",
    "lambda_IcMt_r_info", "lambda_IcMt_nr_info",
    "lambda_ItMt_r_info", "lambda_ItMt_nr_info",
    "last_loglikelihood",
    "n_IcMc", "n_ItMc", "n_IcMt", "n_ItMt",
    "t_judge"
  )
  expect_equal(colnames(estimators), expected_colnames)
})
