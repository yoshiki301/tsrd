calc_induction_rmst <- function (
  estimators,
  t_judge
) {
  # calculate each point estimate of RMST
  rmst_Ic <- (
    estimators$theta_Ic * t_judge +
    (1.0 - stats::dexp(t_judge, estimators$lambda_Ic_nr)) / estimators$lambda_Ic_nr
  )
  rmst_It <- (
    estimators$theta_It * t_judge +
    (1.0 - stats::dexp(t_judge, estimators$lambda_It_nr)) / estimators$lambda_It_nr
  )

  # TODO: calculate CI of each point estimate of RMST
  rmst_Ic_lower_ci <- 0.0
  rmst_Ic_upper_ci <- 1000
  rmst_It_lower_ci <- 0.0
  rmst_It_upper_ci <- 1000

  result <- data.frame(
    rmst_Ic = rmst_Ic,
    rmst_Ic_lower_ci = rmst_Ic_lower_ci,
    rmst_Ic_upper_ci = rmst_Ic_upper_ci,
    rmst_It = rmst_It,
    rmst_It_lower_ci = rmst_It_lower_ci,
    rmst_It_upper_ci = rmst_It_upper_ci
  )
  return (result)
}

calc_maintenance_hr <- function (
    estimators
) {
  # calculate each point estimate of hazard ratio
  hr_Ic <- estimators$lambda_r_IcMt / estimators$lambda_r_IcMc
  hr_It <- estimators$lambda_r_ItMt / estimators$lambda_r_ItMc

  # TODO: calculate CI of each point estimate of hazard ratio
  hr_Ic_lower_ci <- 0.0
  hr_Ic_upper_ci <- 1000
  hr_It_lower_ci <- 0.0
  hr_It_upper_ci <- 1000

  result <- data.frame(
    hr_Ic = hr_Ic,
    hr_Ic_lower_ci = hr_Ic_lower_ci,
    hr_Ic_upper_ci = hr_Ic_upper_ci,
    hr_It = hr_It,
    hr_It_lower_ci = hr_It_lower_ci,
    hr_It_upper_ci = hr_It_upper_ci
  )
  return (result)
}
