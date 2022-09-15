% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimateEM}
\alias{estimateEM}
\title{Parameter estimation using EM algorithm}
\usage{
estimateEM(
  dataset,
  theta_Ic_init = 0.65,
  theta_It_init = 0.75,
  lambda_Ic_nr_init = 0.1,
  lambda_It_nr_init = 0.05,
  lambda_IcMc_r_init = 0.1,
  lambda_IcMc_nr_init = 0.2,
  lambda_ItMc_r_init = 0.08,
  lambda_ItMc_nr_init = 0.16,
  lambda_IcMt_r_init = 0.08,
  lambda_IcMt_nr_init = 0.06,
  lambda_ItMt_r_init = 0.04,
  lambda_ItMt_nr_init = 0.08,
  max_iter = 5000L,
  eps = 1e-05,
  option = NULL
)
}
\description{
This function is to estimate parameters using EM algorithm.
One dataset from generate_scenario is used in estimation.
}