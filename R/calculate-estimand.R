rmst <- function (
  theta,
  lambda,
  t_judge
) {
  return (
    theta * t_judge +
    (1 - theta) * (1.0 - stats::dexp(t_judge, lambda)) / lambda
  )
}

rmst_p_theta <- function (
  theta,
  lambda,
  t_judge
) {
  return (t_judge - (1.0 - stats::dexp(t_judge, lambda)) / lambda)
}

rmst_p_lambda <- function (
  theta,
  lambda,
  t_judge
) {
  return (
    (1 - theta) *
    (stats::dexp(t_judge, lambda) - ((1.0 - stats::dexp(t_judge, lambda)) / lambda^2))
  )
}

calc_induction_rmst <- function (
  estimators,
  dataset,
  alpha = 0.05
) {
  # get fixed data
  n_Ic <- estimators$n_IcMc + estimators$n_IcMt
  n_It <- estimators$n_ItMc + estimators$n_ItMt
  t_judge <- estimators$t_judge

  # get estimated paramters
  theta_Ic <- estimators$theta_Ic
  theta_It <- estimators$theta_It
  lambda_Ic_nr <- estimators$lambda_Ic_nr
  lambda_It_nr <- estimators$lambda_It_nr

  # get variance
  theta_Ic_var <- 1.0 / estimators$theta_Ic_info
  theta_It_var <- 1.0 / estimators$theta_It_info
  lambda_Ic_nr_var <- 1.0 / estimators$lambda_Ic_nr_info
  lambda_It_nr_var <- 1.0 / estimators$lambda_It_nr_info

  # calculate each point estimate of RMST
  rmst_Ic <- rmst(theta_Ic, lambda_Ic_nr, t_judge)
  rmst_It <- rmst(theta_It, lambda_It_nr, t_judge)

  # calculate variance of RMST by delta method
  rmst_Ic_var <- (
    theta_Ic_var * (rmst_p_theta(theta_Ic, lambda_Ic_nr, t_judge))^2 +
    lambda_Ic_nr_var * (rmst_p_lambda(theta_Ic, lambda_Ic_nr, t_judge))^2
  )
  rmst_It_var <- (
    theta_It_var * (rmst_p_theta(theta_It, lambda_It_nr, t_judge))^2 +
    lambda_It_nr_var * (rmst_p_lambda(theta_It, lambda_It_nr, t_judge))^2
  )

  # calculate CI based on asymptotic normality
  rmst_Ic_lower_ci <- rmst_Ic + stats::qnorm(alpha / 2) * sqrt(rmst_Ic_var / n_Ic)
  rmst_Ic_upper_ci <- rmst_Ic + stats::qnorm(1 - (alpha / 2)) * sqrt(rmst_Ic_var / n_Ic)
  rmst_It_lower_ci <- rmst_It + stats::qnorm(alpha / 2) * sqrt(rmst_It_var / n_It)
  rmst_It_upper_ci <- rmst_It + stats::qnorm(1 - (alpha / 2)) * sqrt(rmst_It_var / n_It)

  # calculate true RMST
  theta_Ic_true <- subset(dataset, induction == "Ic", select = theta)[1, 1]
  theta_It_true <- subset(dataset, induction == "It", select = theta)[1, 1]
  lambda_Ic_nr_true <- subset(dataset, induction == "Ic", select = lambda_I_nr)[1, 1]
  lambda_It_nr_true <- subset(dataset, induction == "It", select = lambda_I_nr)[1, 1]
  rmst_Ic_true <- rmst(theta_Ic_true, lambda_Ic_nr_true, t_judge)
  rmst_It_true <- rmst(theta_It_true, lambda_It_nr_true, t_judge)

  result <- data.frame(
    rmst_Ic = rmst_Ic,
    rmst_Ic_lower_ci = rmst_Ic_lower_ci,
    rmst_Ic_upper_ci = rmst_Ic_upper_ci,
    rmst_Ic_true = rmst_Ic_true,
    rmst_It = rmst_It,
    rmst_It_lower_ci = rmst_It_lower_ci,
    rmst_It_upper_ci = rmst_It_upper_ci,
    rmst_It_true = rmst_It_true,
    alpha = alpha
  )
  return (result)
}

hr <- function (
  lambda_IMt,
  lambda_IMc
) {
  return (lambda_IMt / lambda_IMc)
}

hr_p_lambda_IMt <- function (
  lambda_IMt,
  lambda_IMc
) {
  return (1.0 / lambda_IMc)
}

hr_p_lambda_IMc <- function (
  lambda_IMt,
  lambda_IMc
) {
  return (-lambda_IMt / lambda_IMc^2)
}

calc_maintenance_hr <- function (
  estimators,
  dataset,
  alpha = 0.05
) {
  # get fixed data
  n_Ic <- estimators$n_IcMc + estimators$n_IcMt
  n_It <- estimators$n_ItMc + estimators$n_ItMt

  # get estimated paramters
  lambda_IcMc_r <- estimators$lambda_IcMc_r
  lambda_ItMc_r <- estimators$lambda_ItMc_r
  lambda_IcMt_r <- estimators$lambda_IcMt_r
  lambda_ItMt_r <- estimators$lambda_ItMt_r

  # get variance
  lambda_IcMc_r_var <- 1.0 / estimators$lambda_IcMc_r_info
  lambda_ItMc_r_var <- 1.0 / estimators$lambda_ItMc_r_info
  lambda_IcMt_r_var <- 1.0 / estimators$lambda_IcMt_r_info
  lambda_ItMt_r_var <- 1.0 / estimators$lambda_ItMt_r_info

  # calculate each point estimate of hazard ratio
  hr_Ic <- hr(lambda_IcMt_r, lambda_IcMc_r)
  hr_It <- hr(lambda_ItMt_r, lambda_ItMc_r)

  # calculate variance of hazard ratio by delta method
  hr_Ic_var <- (
    lambda_IcMc_r * (hr_p_lambda_IMc(lambda_IcMt_r, lambda_IcMc_r))^2 +
    lambda_IcMt_r * (hr_p_lambda_IMt(lambda_IcMt_r, lambda_IcMc_r))^2
  )
  hr_It_var <- (
    lambda_ItMc_r * (hr_p_lambda_IMc(lambda_ItMt_r, lambda_ItMc_r))^2 +
    lambda_ItMt_r * (hr_p_lambda_IMt(lambda_ItMt_r, lambda_ItMc_r))^2
  )

  # calculate CI based on asymptotic normality
  hr_Ic_lower_ci <- hr_Ic + stats::qnorm(alpha / 2) * sqrt(hr_Ic_var / n_Ic)
  hr_Ic_upper_ci <- hr_Ic + stats::qnorm(1 - (alpha / 2)) * sqrt(hr_Ic_var / n_Ic)
  hr_It_lower_ci <- hr_It + stats::qnorm(alpha / 2) * sqrt(hr_It_var / n_It)
  hr_It_upper_ci <- hr_It + stats::qnorm(1 - (alpha / 2)) * sqrt(hr_It_var / n_It)

  # calculate true hazard ratio
  lambda_IcMc_r_true <- subset(dataset, induction == "Ic" & maintenance == "Mc", select = lambda_IM_r)[1, 1]
  lambda_ItMc_r_true <- subset(dataset, induction == "It" & maintenance == "Mc", select = lambda_IM_r)[1, 1]
  lambda_IcMt_r_true <- subset(dataset, induction == "Ic" & maintenance == "Mt", select = lambda_IM_r)[1, 1]
  lambda_ItMt_r_true <- subset(dataset, induction == "It" & maintenance == "Mt", select = lambda_IM_r)[1, 1]
  hr_Ic_true <- hr(lambda_IcMt_r_true, lambda_IcMc_r_true)
  hr_It_true <- hr(lambda_ItMt_r_true, lambda_ItMc_r_true)

  result <- data.frame(
    hr_Ic = hr_Ic,
    hr_Ic_lower_ci = hr_Ic_lower_ci,
    hr_Ic_upper_ci = hr_Ic_upper_ci,
    hr_Ic_true = hr_Ic_true,
    hr_It = hr_It,
    hr_It_lower_ci = hr_It_lower_ci,
    hr_It_upper_ci = hr_It_upper_ci,
    hr_It_true = hr_It_true,
    alpha = alpha
  )
  return (result)
}
