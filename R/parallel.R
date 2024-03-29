#' Parameter estimation by parallel computing
#'
#' This function is to estimate parameters by parallel computing.
#' Full dataset from generate_scenario is used in estimation.
#' Parallel computing is implemented by parallel package functions.
#'
#' @param dataset The datasets of scenarios from function generate_scenario.
#'
#' @param num_cores The number of assigned cores.
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
estimate_parallely <- function (
  dataset,
  num_cores,
  verbose = T,
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
  # WARNING: this function is about to be removed
  .Deprecated("estimate_sequentially")

  # prepare for parallel computing
  exported_args <- methods::formalArgs(estimateEM)[2:16]
  env <- environment()
  max_cores <- parallel::detectCores()
  if (max_cores < num_cores) {
    stop("num_cores should be less than the value of parallel::detectCores().")
  }
  cl <- parallel::makeCluster(num_cores)
  for (arg_name in exported_args) {
    parallel::clusterExport(cl, arg_name, envir = env)
  }
  parallel::clusterExport(cl, "estimateEM", envir = env)
  parallel::clusterExport(cl, "calc_responsibility", envir = env)
  parallel::clusterExport(cl, "update_theta", envir = env)
  parallel::clusterExport(cl, "update_lambda_I", envir = env)
  parallel::clusterExport(cl, "update_lambda_IM", envir = env)
  parallel::clusterExport(cl, "calc_loglikelihood_stage1", envir = env)
  parallel::clusterExport(cl, "calc_loglikelihood_stage2", envir = env)
  parallel::clusterExport(cl, "calc_Fisher_information", envir = env)
  parallel::clusterExport(cl, "verbose", envir = env)

  parallel::clusterEvalQ(cl, library(testthat)) # for test

  # define function to parallely implement
  target_ids <- unique(dataset$id)
  split_dataset <- lapply(
    target_ids,
    function (target_id) {dataset[dataset$id == target_id,]}
  )
  estimate_func <- function (dataset) {
    if (verbose) {
      # TODO: view progress
    }

    # TODO: enable to select bootstrap
    e <- estimateEM(
      dataset,
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
    estimator <- data.frame(
      # point estimates
      theta_Ic = e$theta_Ic,
      theta_It = e$theta_It,
      lambda_Ic_r = 0.0,
      lambda_Ic_nr = e$lambda_Ic_nr,
      lambda_It_r = 0.0,
      lambda_It_nr = e$lambda_It_nr,
      lambda_IcMc_r = e$lambda_IcMc_r,
      lambda_IcMc_nr = e$lambda_IcMc_nr,
      lambda_ItMc_r = e$lambda_ItMc_r,
      lambda_ItMc_nr = e$lambda_ItMc_nr,
      lambda_IcMt_r = e$lambda_IcMt_r,
      lambda_IcMt_nr = e$lambda_IcMt_nr,
      lambda_ItMt_r = e$lambda_ItMt_r,
      lambda_ItMt_nr = e$lambda_ItMt_nr,
      # Fisher information
      theta_Ic_info = e$theta_Ic_info,
      theta_It_info = e$theta_It_info,
      lambda_Ic_nr_info = e$lambda_Ic_nr_info,
      lambda_It_nr_info = e$lambda_It_nr_info,
      lambda_IcMc_r_info = e$lambda_IcMc_r_info,
      lambda_IcMc_nr_info = e$lambda_IcMc_nr_info,
      lambda_ItMc_r_info = e$lambda_ItMc_r_info,
      lambda_ItMc_nr_info = e$lambda_ItMc_nr_info,
      lambda_IcMt_r_info = e$lambda_IcMt_r_info,
      lambda_IcMt_nr_info = e$lambda_IcMt_nr_info,
      lambda_ItMt_r_info = e$lambda_ItMt_r_info,
      lambda_ItMt_nr_info = e$lambda_ItMt_nr_info
    )
    return (estimator)
  }

  # implement
  estimators <- parallel::parLapply(
    cl = cl, X = split_dataset, fun = estimate_func
  )
  parallel::stopCluster(cl)
  result <- Reduce(rbind, estimators)
  return (as.data.frame(result))
}
