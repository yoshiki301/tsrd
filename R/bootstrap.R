resample_dataset <- function (
  dataset
) {
  # get index of splited dataset into each assignment
  X_IcMc <- subset(dataset, induction == "Ic" & maintenance == "Mc")
  X_ItMc <- subset(dataset, induction == "It" & maintenance == "Mc")
  X_IcMt <- subset(dataset, induction == "Ic" & maintenance == "Mt")
  X_ItMt <- subset(dataset, induction == "It" & maintenance == "Mt")
  index_IcMc <- rownames(X_IcMc)
  index_ItMc <- rownames(X_ItMc)
  index_IcMt <- rownames(X_IcMt)
  index_ItMt <- rownames(X_ItMt)

  # resample dataset
  resampled_IcMc <- X_IcMc[sample(index_IcMc, length(index_IcMc), replace = T),]
  resampled_ItMc <- X_ItMc[sample(index_ItMc, length(index_ItMc), replace = T),]
  resampled_IcMt <- X_IcMt[sample(index_IcMt, length(index_IcMt), replace = T),]
  resampled_ItMt <- X_ItMt[sample(index_ItMt, length(index_ItMt), replace = T),]

  resampled_dataset <- rbind(resampled_IcMc, resampled_ItMc, resampled_IcMt, resampled_ItMt)
  rownames(resampled_dataset) <- NULL # reset index
  return (resampled_dataset)
}

#' Calculation CI based on bootstrap sampling
#'
#' This function is to estimate CI of estimands by bootstrap sampling.
#' Bootstrap CIs are calculated based on asymptotic normality, two-sided.
#' One dataset from generate_scenario is used in estimation.
#'
#' @param dataset The dataset of one scenario from function generate_scenario.
#'
#' @param boot_num The number of bootstrap sampling.
#'
#' @param alpha The significance level.
#'
#' @param theta_Ic_init Initial value vectors of proportion of responders in induction control.
#'
#' @param theta_It_init Initial value vectors of proportion of responders in induction treatment.
#'
#' @param lambda_Ic_nr_init Initial value vectors of hazard of non-responders in induction control.
#'
#' @param lambda_It_nr_init Initial value vectors of hazard of non-responders in induction treatment.
#'
#' @param lambda_IcMc_r_init Initial value vectors of hazard of responders in maintenance control when induction control.
#'
#' @param lambda_IcMc_nr_init Initial value vectors of hazard of non-responders in maintenance control when induction control.
#'
#' @param lambda_ItMc_r_init Initial value vectors of hazard of responders in maintenance treatment when induction control.
#'
#' @param lambda_ItMc_nr_init Initial value vectors of hazard of non-responders in maintenance treatment when induction control.
#'
#' @param lambda_IcMt_r_init Initial value vectors of hazard of responders in maintenance control when induction treatment.
#'
#' @param lambda_IcMt_nr_init Initial value vectors of hazard of non-responders in maintenance control when induction treatment.
#'
#' @param lambda_ItMt_r_init Initial value vectors of hazard of responders in maintenance treatment when induction treatment.
#'
#' @param lambda_ItMt_nr_init Initial value vectors of hazard of non-responders in maintenance treatment when induction treatment.
#'
#' @param max_iter The maximum number of iteration in each EM algorithm.
#'
#' @param eps The threshold of difference log-likelihoods to stop EM algorithm.
#'
#' @param option If NULL, implement default EM algorithm.
#' If "fix_theta", implement EM algorithm being theta fixed.
#'
#' @export
estimateEM_bootstrap_ci <- function (
  dataset,
  boot_num = 1000L,
  alpha = 0.05,
  theta_Ic_init = 0.65,
  theta_It_init = 0.75,
  lambda_Ic_nr_init = 0.10,
  lambda_It_nr_init = 0.05,
  lambda_IcMc_r_init = 0.10,
  lambda_IcMc_nr_init = 0.20,
  lambda_ItMc_r_init = 0.08,
  lambda_ItMc_nr_init = 0.16,
  lambda_IcMt_r_init = 0.08,
  lambda_IcMt_nr_init = 0.06,
  lambda_ItMt_r_init = 0.04,
  lambda_ItMt_nr_init = 0.08,
  max_iter = 5000L,
  eps = 1e-5,
  option = NULL
) {
  boot_func <- function (
    boot_id
  ) {
    resampled_data <- resample_dataset(dataset)
    estimator <- estimateEM_as_frame(
      resampled_data,
      theta_Ic_init = theta_Ic_init,
      theta_It_init = theta_It_init,
      lambda_Ic_nr_init = lambda_Ic_nr_init,
      lambda_It_nr_init = lambda_It_nr_init,
      lambda_IcMc_r_init = lambda_IcMc_r_init,
      lambda_IcMc_nr_init = lambda_IcMc_nr_init,
      lambda_ItMc_r_init = lambda_ItMc_r_init,
      lambda_ItMc_nr_init = lambda_ItMc_nr_init,
      lambda_IcMt_r_init = lambda_IcMt_r_init,
      lambda_IcMt_nr_init = lambda_IcMt_nr_init,
      lambda_ItMt_r_init = lambda_ItMt_r_init,
      lambda_ItMt_nr_init = lambda_ItMt_nr_init,
      max_iter = max_iter,
      eps = eps,
      option = option
    )

    rmst <- calc_induction_rmst(estimator, resampled_data, alpha = alpha)
    hr <- calc_maintenance_hr(estimator, resampled_data, alpha = alpha)

    estimands <- data.frame(
      rmst_Ic = rmst$rmst_Ic,
      rmst_Ic_true = rmst$rmst_Ic_true,
      rmst_It = rmst$rmst_It,
      rmst_It_true = rmst$rmst_It_true,
      hr_Ic = hr$hr_Ic,
      hr_Ic_true = hr$hr_Ic_true,
      hr_It = hr$hr_It,
      hr_It_true = hr$hr_It_true
    )
    return (estimands)
  }

  # implement bootstrap
  boot_ids <- 1:boot_num
  boot_list <- lapply(boot_ids, boot_func)
  boot_estimators <- Reduce(rbind, boot_list)

  # two-sided bootstrap CI
  lower_percentile <- alpha / 2
  upper_percentile <- 1 - (alpha / 2)
  lower_quantile <- sapply(boot_estimators, function (x) stats::quantile(x, lower_percentile))
  upper_quantile <- sapply(boot_estimators, function (x) stats::quantile(x, upper_percentile))

  percentile <- c(lower_percentile, upper_percentile)
  result <- cbind(rbind(lower_quantile, upper_quantile), percentile)
  colnames(result) <- c(
    "rmst_Ic", "rmst_Ic_true", "rmst_It", "rmst_It_true",
    "hr_Ic", "hr_Ic_true", "hr_It", "hr_It_true",
    "percentile"
  )
  rownames(result) <- NULL
  return (as.data.frame(result))
}

#' Calculation CI by sequential execution
#'
#' This function is to estimate CI of estimands by sequential execution.
#' Full dataset from generate_scenario is used in estimation.
#'
#' @param dataset The datasets of scenarios from function generate_scenario.
#'
#' @param boot_num The number of bootstrap sampling.
#'
#' @param alpha The significance level.
#'
#' @param verbose If TRUE, show progress.
#'
#' @param theta_Ic_init Initial value vectors of proportion of responders in induction control.
#'
#' @param theta_It_init Initial value vectors of proportion of responders in induction treatment.
#'
#' @param lambda_Ic_nr_init Initial value vectors of hazard of non-responders in induction control.
#'
#' @param lambda_It_nr_init Initial value vectors of hazard of non-responders in induction treatment.
#'
#' @param lambda_IcMc_r_init Initial value vectors of hazard of responders in maintenance control when induction control.
#'
#' @param lambda_IcMc_nr_init Initial value vectors of hazard of non-responders in maintenance control when induction control.
#'
#' @param lambda_ItMc_r_init Initial value vectors of hazard of responders in maintenance treatment when induction control.
#'
#' @param lambda_ItMc_nr_init Initial value vectors of hazard of non-responders in maintenance treatment when induction control.
#'
#' @param lambda_IcMt_r_init Initial value vectors of hazard of responders in maintenance control when induction treatment.
#'
#' @param lambda_IcMt_nr_init Initial value vectors of hazard of non-responders in maintenance control when induction treatment.
#'
#' @param lambda_ItMt_r_init Initial value vectors of hazard of responders in maintenance treatment when induction treatment.
#'
#' @param lambda_ItMt_nr_init Initial value vectors of hazard of non-responders in maintenance treatment when induction treatment.
#'
#' @param max_iter The maximum number of iteration in each EM algorithm.
#'
#' @param eps The threshold of difference log-likelihoods to stop EM algorithm.
#'
#' @param option If NULL, implement default EM algorithm.
#' If "fix_theta", implement EM algorithm being theta fixed.
#'
#' @export
estimate_boot_ci_sequentially <- function (
  dataset,
  boot_num = 1000L,
  alpha = 0.05,
  verbose = TRUE,
  theta_Ic_init = 0.65,
  theta_It_init = 0.75,
  lambda_Ic_nr_init = 0.10,
  lambda_It_nr_init = 0.05,
  lambda_IcMc_r_init = 0.10,
  lambda_IcMc_nr_init = 0.20,
  lambda_ItMc_r_init = 0.08,
  lambda_ItMc_nr_init = 0.16,
  lambda_IcMt_r_init = 0.08,
  lambda_IcMt_nr_init = 0.06,
  lambda_ItMt_r_init = 0.04,
  lambda_ItMt_nr_init = 0.08,
  max_iter = 5000L,
  eps = 1e-5,
  option = NULL
) {
  target_ids <- unique(dataset$id)
  split_dataset <- lapply(
    target_ids,
    function (target_id) {dataset[dataset$id == target_id,]}
  )
  ci_func <- function (dataset) {
    target_id <- dataset$id[1]
    if (verbose) {
      cat(target_id, fill = TRUE)
    }
    confidence_interval <- estimateEM_bootstrap_ci(
      dataset,
      boot_num = boot_num,
      alpha = alpha,
      theta_Ic_init = theta_Ic_init,
      theta_It_init = theta_Ic_init,
      lambda_Ic_nr_init = lambda_Ic_nr_init,
      lambda_It_nr_init = lambda_It_nr_init,
      lambda_IcMc_r_init = lambda_IcMc_r_init,
      lambda_IcMc_nr_init = lambda_IcMc_nr_init,
      lambda_ItMc_r_init = lambda_ItMc_r_init,
      lambda_ItMc_nr_init = lambda_ItMc_nr_init,
      lambda_IcMt_r_init = lambda_IcMt_r_init,
      lambda_IcMt_nr_init = lambda_IcMt_nr_init,
      lambda_ItMt_r_init = lambda_ItMt_r_init,
      lambda_ItMt_nr_init = lambda_ItMt_nr_init,
      max_iter = max_iter,
      eps = eps,
      option = option
    )
    cbind(
      id = target_id,
      confidence_interval
    )
  }

  # execution
  confidence_intervals <- lapply(
    split_dataset, ci_func
  )
  result <- Reduce(rbind, confidence_intervals)
  result
}

#' Measure the coverage probability
#'
#' This function is to measure the coverage probability.
#'
#' @param ci The result of CI calculation.
#'
#' @export
coverage <- function (
  ci
) {
  target_ids <- unique(ci$id)
  upper <- max(ci$percentile)
  lower <- min(ci$percentile)
  coverage_func <- function (target_id) {
    ci_upper <- subset(ci, id == target_id & percentile == upper)
    ci_lower <- subset(ci, id == target_id & percentile == lower)

    rmst_Ic_covered <- ifelse(
      ci_lower$rmst_Ic <= ci_lower$rmst_Ic_true & ci_upper$rmst_Ic_true <= ci_upper$rmst_Ic,
      1, 0
    )
    rmst_Ic_lower <- ifelse(ci_lower$rmst_Ic > ci_lower$rmst_Ic_true, 1, 0)
    rmst_Ic_upper <- ifelse(ci_upper$rmst_Ic < ci_upper$rmst_Ic_true, 1, 0)
    rmst_It_covered <- ifelse(
      ci_lower$rmst_It <= ci_lower$rmst_It_true & ci_upper$rmst_It_true <= ci_upper$rmst_It,
      1, 0
    )
    rmst_It_lower <- ifelse(ci_lower$rmst_It > ci_lower$rmst_It_true, 1, 0)
    rmst_It_upper <- ifelse(ci_upper$rmst_It < ci_upper$rmst_It_true, 1, 0)
    hr_Ic_covered <- ifelse(
      ci_lower$hr_Ic <= ci_lower$hr_Ic_true & ci_upper$hr_Ic_true <= ci_upper$hr_Ic,
      1, 0
    )
    hr_Ic_lower <- ifelse(ci_lower$hr_Ic > ci_lower$hr_Ic_true, 1, 0)
    hr_Ic_upper <- ifelse(ci_upper$hr_Ic < ci_upper$hr_Ic_true, 1, 0)
    hr_It_covered <- ifelse(
      ci_lower$hr_It <= ci_lower$hr_It_true & ci_upper$hr_It_true <= ci_upper$hr_It,
      1, 0
    )
    hr_It_lower <- ifelse(ci_lower$hr_It > ci_lower$hr_It_true, 1, 0)
    hr_It_upper <- ifelse(ci_upper$hr_It < ci_upper$hr_It_true, 1, 0)

    data.frame(
      id = target_id,
      rmst_Ic_covered = rmst_Ic_covered,
      rmst_Ic_lower = rmst_Ic_lower,
      rmst_Ic_upper = rmst_Ic_upper,
      rmst_It_covered = rmst_It_covered,
      rmst_It_lower = rmst_It_lower,
      rmst_It_upper = rmst_It_upper,
      hr_Ic_covered = hr_Ic_covered,
      hr_Ic_lower = hr_Ic_lower,
      hr_Ic_upper = hr_Ic_upper,
      hr_It_covered = hr_It_covered,
      hr_It_lower = hr_It_lower,
      hr_It_upper = hr_It_upper
    )
  }

  coverage_list <- lapply(target_ids, coverage_func)
  Reduce(rbind, coverage_list)
}
