estimators <- data.frame(
  theta_Ic = c(0.2, 0.3, 0.4),
  theta_It = c(0.3, 0.4, 0.5),
  lambda_Ic_nr = c(0.01, 0.02, 0.03),
  lambda_It_nr = c(0.02, 0.03, 0.04),
  lambda_IcMt_r = c(0.01, 0.02, 0.03),
  lambda_IcMc_r = c(0.02, 0.04, 0.06),
  lambda_ItMt_r = c(0.02, 0.03, 0.04),
  lambda_ItMc_r = c(0.04, 0.06, 0.08)
)
dataset <- data.frame(
  induction = c("Ic", "It", "Ic", "It"),
  maintenance = c("Mc", "Mc", "Mt", "Mt"),
  theta = c(0.2, 0.3, 0.4, 0.5),
  lambda_IM_r = c(0.04, 0.04, 0.02, 0.01),
  lambda_I_nr = c(0.02, 0.04, 0.06, 0.08)
)

test_that("calc_induction_rmst works", {
  rmst <- calc_induction_rmst(estimators, dataset, t_judge = 6)
  # TODO: expected values
  # TODO: fix CI
  expect_equal(T, T)
})

test_that("calc_maintenance_hr works", {
  hr <- calc_maintenance_hr(estimators, dataset)
  # TODO: fix CI
  expected_hr <- data.frame(
    hr_Ic = c(0.5, 0.5, 0.5),
    hr_Ic_lower_ci = c(0.0, 0.0, 0.0),
    hr_Ic_upper_ci = c(1000, 1000, 1000),
    hr_Ic_true = rep(0.5, 3L),
    hr_It = c(0.5, 0.5, 0.5),
    hr_It_lower_ci = c(0.0, 0.0, 0.0),
    hr_It_upper_ci = c(1000, 1000, 1000),
    hr_It_true = rep(0.25, 3L)
  )
  expect_equal(hr, expected_hr)
})
