rmst <- function (
  theta,
  lambda,
  t_judge
) {
  theta * t_judge + (1 - theta) * (1.0 - stats::dexp(t_judge, lambda)) / lambda
}

rmst_p_theta <- function (
  theta,
  lambda,
  t_judge
) {
  t_judge - (1.0 - stats::dexp(t_judge, lambda)) / lambda
}

rmst_p_lambda <- function (
  theta,
  lambda,
  t_judge
) {
    (1 - theta) * (stats::dexp(t_judge, lambda) - ((1.0 - stats::dexp(t_judge, lambda)) / lambda^2))
}

rmst_p2_theta_lambda <- function (
  theta,
  lambda,
  t_judge
) {
  -(stats::dexp(t_judge, lambda) - ((1.0 - stats::dexp(t_judge, lambda)) / lambda^2))
}

calc_induction_rmst <- function (
  estimators,
  dataset,
  alpha = 0.05,
  monte_carlo_num = 100L,
  pseudo_data_censor_rate = 0.0,
  seed = 42L
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

  # get variance matrix via Louis' method
  pseudo_datasets <- generate_pseudo_dataset(
    original_dataset = dataset,
    estimators = estimators,
    num = monte_carlo_num,
    censor_rate = pseudo_data_censor_rate,
    seed = seed
  )
  fisher_infomation <- fisher_information_monte_carlo(pseudo_datasets, estimators)
  variance_matrix <- solve(fisher_infomation, tol = 1e-100)

  # get variance/covariance associated with RMST
  theta_Ic_var <- variance_matrix[1,1]
  theta_It_var <- variance_matrix[2,2]
  lambda_Ic_nr_var <- variance_matrix[3,3]
  lambda_It_nr_var <- variance_matrix[4,4]
  theta_Ic_lambda_Ic_nr_cov <- variance_matrix[1,3]
  theta_It_lambda_It_nr_cov <- variance_matrix[2,4]
  #theta_Ic_lambda_Ic_nr_cov <- 0.0
  #theta_It_lambda_It_nr_cov <- 0.0

  # calculate each point estimate of RMST
  rmst_Ic <- rmst(theta_Ic, lambda_Ic_nr, t_judge)
  rmst_It <- rmst(theta_It, lambda_It_nr, t_judge)

  # calculate variance of RMST by delta method
  rmst_Ic_var <- abs(
    theta_Ic_var * (rmst_p_theta(theta_Ic, lambda_Ic_nr, t_judge))^2 +
    lambda_Ic_nr_var * (rmst_p_lambda(theta_Ic, lambda_Ic_nr, t_judge))^2 +
    2 * theta_Ic_lambda_Ic_nr_cov * rmst_p2_theta_lambda(theta_Ic, lambda_Ic_nr, t_judge)
  )
  rmst_It_var <- abs(
    theta_It_var * (rmst_p_theta(theta_It, lambda_It_nr, t_judge))^2 +
    lambda_It_nr_var * (rmst_p_lambda(theta_It, lambda_It_nr, t_judge))^2 +
    2 * theta_It_lambda_It_nr_cov * rmst_p2_theta_lambda(theta_It, lambda_It_nr, t_judge)
  )

  # calculate CI based on asymptotic normality
  rmst_Ic_lower_ci <- rmst_Ic + stats::qnorm(alpha / 2) * sqrt(rmst_Ic_var)
  rmst_Ic_upper_ci <- rmst_Ic + stats::qnorm(1 - (alpha / 2)) * sqrt(rmst_Ic_var)
  rmst_It_lower_ci <- rmst_It + stats::qnorm(alpha / 2) * sqrt(rmst_It_var)
  rmst_It_upper_ci <- rmst_It + stats::qnorm(1 - (alpha / 2)) * sqrt(rmst_It_var)

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
  lambda_IMt / lambda_IMc
}

hr_p_lambda_IMt <- function (
  lambda_IMt,
  lambda_IMc
) {
  1.0 / lambda_IMc
}

hr_p_lambda_IMc <- function (
  lambda_IMt,
  lambda_IMc
) {
  -lambda_IMt / lambda_IMc^2
}

hr_p2_lambda_IMt_lambda_IMc <- function (
  lambda_IMt,
  lambda_IMc
) {
  -1.0 / lambda_IMc^2
}

calc_maintenance_hr <- function (
  estimators,
  dataset,
  alpha = 0.05,
  monte_carlo_num = 100L,
  pseudo_data_censor_rate = 0.0,
  seed = 42L
) {
  # get fixed data
  n_Ic <- estimators$n_IcMc + estimators$n_IcMt
  n_It <- estimators$n_ItMc + estimators$n_ItMt

  # get estimated paramters
  lambda_IcMc_r <- estimators$lambda_IcMc_r
  lambda_ItMc_r <- estimators$lambda_ItMc_r
  lambda_IcMt_r <- estimators$lambda_IcMt_r
  lambda_ItMt_r <- estimators$lambda_ItMt_r

  # get variance matrix via Louis' method
  pseudo_datasets <- generate_pseudo_dataset(
    original_dataset = dataset,
    estimators = estimators,
    num = monte_carlo_num,
    censor_rate = pseudo_data_censor_rate,
    seed = seed
  )
  fisher_infomation <- fisher_information_monte_carlo(pseudo_datasets, estimators)
  variance_matrix <- solve(fisher_infomation, tol = 1e-100)

  # get variance/covariance associated with RMST
  lambda_IcMc_r_var <- variance_matrix[5,5]
  lambda_ItMc_r_var <- variance_matrix[7,7]
  lambda_IcMt_r_var <- variance_matrix[9,9]
  lambda_ItMt_r_var <- variance_matrix[11,11]
  lambda_IcMc_r_lambda_IcMt_r_cov <- variance_matrix[5,9]
  lambda_ItMc_r_lambda_ItMt_r_cov <- variance_matrix[7,11]
  #lambda_IcMc_r_lambda_IcMt_r_cov <- 0.0
  #lambda_ItMc_r_lambda_ItMt_r_cov <- 0.0

  # calculate each point estimate of hazard ratio
  hr_Ic <- hr(lambda_IcMt_r, lambda_IcMc_r)
  hr_It <- hr(lambda_ItMt_r, lambda_ItMc_r)

  # calculate variance of hazard ratio by delta method
  hr_Ic_var <- abs(
    lambda_IcMc_r_var * (hr_p_lambda_IMc(lambda_IcMt_r, lambda_IcMc_r))^2 +
    lambda_IcMt_r_var * (hr_p_lambda_IMt(lambda_IcMt_r, lambda_IcMc_r))^2 +
    2 * lambda_IcMc_r_lambda_IcMt_r_cov * hr_p2_lambda_IMt_lambda_IMc(lambda_IcMt_r, lambda_IcMc_r)
  )
  hr_It_var <- abs(
    lambda_ItMc_r_var * (hr_p_lambda_IMc(lambda_ItMt_r, lambda_ItMc_r))^2 +
    lambda_ItMt_r_var * (hr_p_lambda_IMt(lambda_ItMt_r, lambda_ItMc_r))^2 +
    2 * lambda_ItMc_r_lambda_ItMt_r_cov * hr_p2_lambda_IMt_lambda_IMc(lambda_ItMt_r, lambda_ItMc_r)
  )

  # calculate CI based on asymptotic normality
  hr_Ic_lower_ci <- hr_Ic + stats::qnorm(alpha / 2) * sqrt(hr_Ic_var)
  hr_Ic_upper_ci <- hr_Ic + stats::qnorm(1 - (alpha / 2)) * sqrt(hr_Ic_var)
  hr_It_lower_ci <- hr_It + stats::qnorm(alpha / 2) * sqrt(hr_It_var)
  hr_It_upper_ci <- hr_It + stats::qnorm(1 - (alpha / 2)) * sqrt(hr_It_var)

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

calc_treatment_effects_sequentially <- function (
  estimators,
  datasets,
  verbose = TRUE,
  alpha = 0.05,
  monte_carlo_num = 100L,
  pseudo_data_censor_rate = 0.0,
  seed = 42L
) {
  ids <- unique(estimators$id)
  effect_func <- function (target_id) {
    if (verbose) {
      cat(target_id, fill = TRUE)
    }
    estimator <- estimators[estimators$id == target_id,]
    dataset <- datasets[datasets$id == target_id,]
    rmst <- calc_induction_rmst(
      estimators = estimator, dataset = dataset, alpha = alpha,
      monte_carlo_num = monte_carlo_num, pseudo_data_censor_rate = pseudo_data_censor_rate,
      seed = seed
    )
    hr <- calc_maintenance_hr(
      estimators = estimator, dataset = dataset, alpha = alpha,
      monte_carlo_num = monte_carlo_num, pseudo_data_censor_rate = pseudo_data_censor_rate,
      seed = seed
    )
    data.frame(
      id = target_id,
      rmst_Ic = rmst$rmst_Ic,
      rmst_Ic_lower = rmst$rmst_Ic_lower_ci,
      rmst_Ic_upper = rmst$rmst_Ic_upper_ci,
      rmst_Ic_true = rmst$rmst_Ic_true,
      rmst_It = rmst$rmst_It,
      rmst_It_lower = rmst$rmst_It_lower_ci,
      rmst_It_upper = rmst$rmst_It_upper_ci,
      rmst_It_true = rmst$rmst_It_true,
      hr_Ic = hr$hr_Ic,
      hr_Ic_lower = hr$hr_Ic_lower_ci,
      hr_Ic_upper = hr$hr_Ic_upper_ci,
      hr_Ic_true = hr$hr_Ic_true,
      hr_It = hr$hr_It,
      hr_It_lower = hr$hr_It_lower_ci,
      hr_It_upper = hr$hr_It_upper_ci,
      hr_It_true = hr$hr_It_true,
      alpha = alpha
    )
  }

  effect_list <- lapply(ids, effect_func)
  Reduce(rbind, effect_list)
}

calc_bias <- function (
  treatment_effects
) {
  data.frame(
    rmst_Ic_bias = treatment_effects$rmst_Ic - treatment_effects$rmst_Ic_true,
    rmst_It_bias = treatment_effects$rmst_It - treatment_effects$rmst_It_true,
    hr_Ic_bias = treatment_effects$hr_Ic - treatment_effects$hr_Ic_true,
    hr_It_bias = treatment_effects$hr_It - treatment_effects$hr_It_true
  )
}

calc_coverage_prob <- function (
  treatment_effects
) {
  coverage_frame <- data.frame(
    rmst_Ic_covered = ifelse(
      (treatment_effects$rmst_Ic_lower <= treatment_effects$rmst_Ic_true
       & treatment_effects$rmst_Ic_upper >= treatment_effects$rmst_Ic_true), 1, 0
    ),
    rmst_Ic_under = ifelse(
      treatment_effects$rmst_Ic_lower > treatment_effects$rmst_Ic_true, 1, 0
    ),
    rmst_Ic_over = ifelse(
      treatment_effects$rmst_Ic_upper < treatment_effects$rmst_Ic_true, 1, 0
    ),
    rmst_It_covered = ifelse(
      (treatment_effects$rmst_It_lower <= treatment_effects$rmst_It_true
       & treatment_effects$rmst_It_upper >= treatment_effects$rmst_It_true), 1, 0
    ),
    rmst_It_under = ifelse(
      treatment_effects$rmst_It_lower > treatment_effects$rmst_It_true, 1, 0
    ),
    rmst_It_over = ifelse(
      treatment_effects$rmst_It_upper < treatment_effects$rmst_It_true, 1, 0
    ),
    hr_Ic_covered = ifelse(
      (treatment_effects$hr_Ic_lower <= treatment_effects$hr_Ic_true
       & treatment_effects$hr_Ic_upper >= treatment_effects$hr_Ic_true), 1, 0
    ),
    hr_Ic_under = ifelse(
      treatment_effects$hr_Ic_lower > treatment_effects$hr_Ic_true, 1, 0
    ),
    hr_Ic_over = ifelse(
      treatment_effects$hr_Ic_upper < treatment_effects$hr_Ic_true, 1, 0
    ),
    hr_It_covered = ifelse(
      (treatment_effects$hr_It_lower <= treatment_effects$hr_It_true
       & treatment_effects$hr_It_upper >= treatment_effects$hr_It_true), 1, 0
    ),
    hr_It_under = ifelse(
      treatment_effects$hr_It_lower > treatment_effects$hr_It_true, 1, 0
    ),
    hr_It_over = ifelse(
      treatment_effects$hr_It_upper < treatment_effects$hr_It_true, 1, 0
    )
  )
  colMeans(coverage_frame)
}
