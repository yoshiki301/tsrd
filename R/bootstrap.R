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
    estimator <- estimateEM(
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

    rmst <- calc_induction_rmst(estimator, dataset, dataset$t_judge[1])
    hr <- calc_maintenance_hr(estimator, dataset)

    estimands <- data.frame(
      rmst_Ic = rmst$rmst_Ic,
      rmst_It = rmst$rmst_It,
      hr_Ic = hr$hr_Ic,
      hr_It = hr$hr_It
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
  colnames(result) <- c("rmst_Ic", "rmst_It", "hr_Ic", "hr_It", "percentile")
  rownames(result) <- NULL
  return (as.data.frame(result))
}
