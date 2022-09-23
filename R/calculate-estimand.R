calc_induction_rmst <- function (
  estimators,
  dataset,
  t_judge
) {
  # calculate each point estimate of RMST
  rmst_Ic <- (
    estimators$theta_Ic * t_judge +
    (1.0 - estimators$theta_Ic) * (1.0 - stats::dexp(t_judge, estimators$lambda_Ic_nr)) / estimators$lambda_Ic_nr
  )
  rmst_It <- (
    estimators$theta_It * t_judge +
    (1.0 - estimators$theta_It) * (1.0 - stats::dexp(t_judge, estimators$lambda_It_nr)) / estimators$lambda_It_nr
  )

  # TODO: calculate CI of each point estimate of RMST
  rmst_Ic_lower_ci <- 0.0
  rmst_Ic_upper_ci <- 1000
  rmst_It_lower_ci <- 0.0
  rmst_It_upper_ci <- 1000

  # calculate true RMST
  theta_Ic_true <- subset(dataset, induction == "Ic", select = theta)[1, 1]
  theta_It_true <- subset(dataset, induction == "It", select = theta)[1, 1]
  lambda_Ic_nr_true <- subset(dataset, induction == "Ic", select = lambda_I_nr)[1, 1]
  lambda_It_nr_true <- subset(dataset, induction == "It", select = lambda_I_nr)[1, 1]
  rmst_Ic_true <- (
    theta_Ic_true * t_judge +
    (1.0 - theta_Ic_true) * (1.0 - stats::dexp(t_judge, lambda_Ic_nr_true)) / lambda_Ic_nr_true
  )
  rmst_It_true <- (
    theta_It_true * t_judge +
    (1.0 - theta_It_true) * (1.0 - stats::dexp(t_judge, lambda_It_nr_true)) / lambda_It_nr_true
  )

  result <- data.frame(
    rmst_Ic = rmst_Ic,
    rmst_Ic_lower_ci = rmst_Ic_lower_ci,
    rmst_Ic_upper_ci = rmst_Ic_upper_ci,
    rmst_Ic_true = rmst_Ic_true,
    rmst_It = rmst_It,
    rmst_It_lower_ci = rmst_It_lower_ci,
    rmst_It_upper_ci = rmst_It_upper_ci,
    rmst_It_true = rmst_It_true
  )
  return (result)
}

calc_maintenance_hr <- function (
  estimators,
  dataset
) {
  # calculate each point estimate of hazard ratio
  hr_Ic <- estimators$lambda_r_IcMt / estimators$lambda_r_IcMc
  hr_It <- estimators$lambda_r_ItMt / estimators$lambda_r_ItMc

  # TODO: calculate CI of each point estimate of hazard ratio
  hr_Ic_lower_ci <- 0.0
  hr_Ic_upper_ci <- 1000
  hr_It_lower_ci <- 0.0
  hr_It_upper_ci <- 1000

  # calculate true hazard ratio
  lambda_r_IcMc_true <- subset(dataset, induction == "Ic" & maintenance == "Mc", select = lambda_IM_r)[1, 1]
  lambda_r_ItMc_true <- subset(dataset, induction == "It" & maintenance == "Mc", select = lambda_IM_r)[1, 1]
  lambda_r_IcMt_true <- subset(dataset, induction == "Ic" & maintenance == "Mt", select = lambda_IM_r)[1, 1]
  lambda_r_ItMt_true <- subset(dataset, induction == "It" & maintenance == "Mt", select = lambda_IM_r)[1, 1]
  hr_Ic_true <- lambda_r_IcMt_true / lambda_r_IcMc_true
  hr_It_true <- lambda_r_ItMt_true / lambda_r_ItMc_true

  result <- data.frame(
    hr_Ic = hr_Ic,
    hr_Ic_lower_ci = hr_Ic_lower_ci,
    hr_Ic_upper_ci = hr_Ic_upper_ci,
    hr_Ic_true = hr_Ic_true,
    hr_It = hr_It,
    hr_It_lower_ci = hr_It_lower_ci,
    hr_It_upper_ci = hr_It_upper_ci,
    hr_It_true = hr_It_true
  )
  return (result)
}
