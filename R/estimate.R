calc_responsibility <- function (
  X_IM,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  time <- X_IM$time
  event <- X_IM$cens
  t_judge <- X_IM$t_judge
  is_I <- X_IM$is_induction

  S_r_judge <- 1 - stats::pexp(t_judge, 0)
  S_nr_judge <- 1 - stats::pexp(t_judge, lambda_I_nr)
  f_r <- (
    is_I * 0.0 +
    (1 - is_I) * S_r_judge * stats::dexp(time - t_judge, lambda_IM_r)
  )
  f_nr <- (
    is_I * stats::dexp(time, lambda_I_nr) +
    (1 - is_I) * S_nr_judge * stats::dexp(time - t_judge, lambda_IM_nr)
  )
  S_r <- (
    is_I * 1.0 +
    (1 - is_I) * S_r_judge * (1 - stats::pexp(time - t_judge, lambda_IM_r))
  )
  S_nr <- (
    is_I * (1 - stats::pexp(time, lambda_I_nr)) +
    (1 - is_I) * S_nr_judge * (1 - stats::pexp(time - t_judge, lambda_IM_nr))
  )

  observation <- theta * f_r / (
    theta * f_r + (1 - theta) * f_nr
  )
  missing <- theta * S_r / (
    theta * S_r + (1 - theta) * S_nr
  )
  responsibility <- event * observation + (1 - event) * missing
  return (responsibility)
}

update_theta <- function (
  X_IMc,
  X_IMt,
  pi_IMc,
  pi_IMt
) {
  theta <- mean(c(pi_IMc, pi_IMt))
  return (theta)
}

update_lambda_I <- function (
  X_IMc,
  X_IMt,
  pi_IMc,
  pi_IMt
) {
  X_I <- rbind(X_IMc, X_IMt)
  pi_I <- c(pi_IMc, pi_IMt)
  time <- X_I$time
  event <- X_I$cens
  t_judge <- X_I$t_judge
  is_I <- X_I$is_induction
  lambda <- sum((1 - pi_I) * is_I * event) / sum(
    (1 - pi_I) * (is_I * time + (1 - is_I) * t_judge)
  )
  return (lambda)
}

update_lambda_IM <- function (
  X_IM,
  pi_IM
) {
  time <- X_IM$time
  event <- X_IM$cens
  t_judge <- X_IM$t_judge
  is_I <- X_IM$is_induction
  lambda_r <- sum((1 - is_I) * event * pi_IM) / sum((1 - is_I) * pi_IM * (time - t_judge))
  lambda_nr <- sum((1 - is_I) * event * (1 - pi_IM)) / sum((1 - is_I) * (1 - pi_IM) * (time - t_judge))
  return (list(
    r = lambda_r,
    nr = lambda_nr
  ))
}

calc_loglikelihood_stage1 <- function (
  X_IMc,
  X_IMt,
  theta,
  lambda_I_nr
) {
  X_I <- rbind(X_IMc, X_IMt)
  time <- X_I$time
  is_I <- X_I$is_induction
  loglikelihood <- sum(is_I) * log(1 - theta) + sum(is_I * stats::dexp(time, lambda_I_nr))
  return (loglikelihood)
}

calc_loglikelihood_stage2 <- function (
  X_IM,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  time <- X_IM$time
  t_judge <- X_IM$t_judge
  is_I <- X_IM$is_induction
  S_r_judge <- 1 - stats::pexp(t_judge, 0)
  S_nr_judge <- 1 - stats::pexp(t_judge, lambda_I_nr)
  loglikelihood <- sum((1 - is_I) * log(
    theta * S_r_judge * stats::dexp(time - t_judge, lambda_IM_r) +
    (1 - theta) * S_nr_judge * stats::dexp(time - t_judge, lambda_IM_nr)
  ))
  loglikelihood[is.nan(loglikelihood)] <- 0.0 # fill induction's log-likelihoods
  return (loglikelihood)
}

calc_Fisher_information <- function (
  dataset,
  parameters
) {
  # fixed data
  X_IcMc <- subset(
    dataset, induction == "Ic" & maintenance == "Mc",
    select = c(cens, is_induction)
  )
  X_ItMc <- subset(
    dataset, induction == "It" & maintenance == "Mc",
    select = c(cens, is_induction)
  )
  X_IcMt <- subset(
    dataset, induction == "Ic" & maintenance == "Mt",
    select = c(cens, is_induction)
  )
  X_ItMt <- subset(
    dataset, induction == "It" & maintenance == "Mt",
    select = c(cens, is_induction)
  )

  X_Ic <- rbind(X_IcMc, X_IcMt)
  X_It <- rbind(X_ItMc, X_ItMt)
  pi_Ic <- c(parameters$pi_IcMc, parameters$pi_IcMt)
  pi_It <- c(parameters$pi_ItMc, parameters$pi_ItMt)

  # calculate the Fisher information of theta
  theta_Ic_info <- sum(
    (pi_Ic / parameters$theta_Ic^2) +
    ((1 - pi_Ic) / (1 - parameters$theta_Ic)^2)
  )
  theta_It_info <- sum(
    (pi_It / parameters$theta_It^2) +
    ((1 - pi_It) / (1 - parameters$theta_It)^2)
  )

  # calculate the Fisher information of hazard of non-responders in induction
  lambda_Ic_nr_info <- sum(
    (1 - pi_Ic) * X_Ic$cens * X_Ic$is_induction / (parameters$lambda_Ic_nr^2)
  )
  lambda_It_nr_info <- sum(
    (1 - pi_It) * X_It$cens * X_It$is_induction / (parameters$lambda_It_nr^2)
  )

  # calculate the Fisher information of hazard in maintenance
  lambda_IcMc_r_info <- sum(
    parameters$pi_IcMc * X_IcMc$cens * (1 - X_IcMc$is_induction) / (parameters$lambda_IcMc_r^2)
  )
  lambda_IcMc_nr_info <- sum(
    (1 - parameters$pi_IcMc) * X_IcMc$cens * (1 - X_IcMc$is_induction) / (parameters$lambda_IcMc_nr^2)
  )
  lambda_ItMc_r_info <- sum(
    parameters$pi_ItMc * X_ItMc$cens * (1 - X_ItMc$is_induction) / (parameters$lambda_ItMc_r^2)
  )
  lambda_ItMc_nr_info <- sum(
    (1 - parameters$pi_ItMc) * X_ItMc$cens * (1 - X_ItMc$is_induction) / (parameters$lambda_ItMc_nr^2)
  )
  lambda_IcMt_r_info <- sum(
    parameters$pi_IcMt * X_IcMt$cens * (1 - X_IcMt$is_induction) / (parameters$lambda_IcMt_r^2)
  )
  lambda_IcMt_nr_info <- sum(
    (1 - parameters$pi_IcMt) * X_IcMt$cens * (1 - X_IcMt$is_induction) / (parameters$lambda_IcMt_nr^2)
  )
  lambda_ItMt_r_info <- sum(
    parameters$pi_ItMt * X_ItMt$cens * (1 - X_ItMt$is_induction) / (parameters$lambda_ItMt_r^2)
  )
  lambda_ItMt_nr_info <- sum(
    (1 - parameters$pi_ItMt) * X_ItMt$cens * (1 - X_ItMt$is_induction) / (parameters$lambda_ItMt_nr^2)
  )

  # TODO: calculate non-diagonal elements

  result <- list(
    theta_Ic_info = theta_Ic_info,
    theta_It_info = theta_It_info,
    lambda_Ic_nr_info = lambda_Ic_nr_info,
    lambda_It_nr_info = lambda_It_nr_info,
    lambda_IcMc_r_info = lambda_IcMc_r_info,
    lambda_IcMc_nr_info = lambda_IcMc_nr_info,
    lambda_ItMc_r_info = lambda_ItMc_r_info,
    lambda_ItMc_nr_info = lambda_ItMc_nr_info,
    lambda_IcMt_r_info = lambda_IcMt_r_info,
    lambda_IcMt_nr_info = lambda_IcMt_nr_info,
    lambda_ItMt_r_info = lambda_ItMt_r_info,
    lambda_ItMt_nr_info = lambda_ItMt_nr_info
  )
  return (result)
}

calc_parameter_ci <- function (
  dataset,
  paramters,
  alpha = 0.05,
  full_denom = TRUE
) {
  fisher_info <- calc_Fisher_information(dataset, parameters)

  # CI by normal approximation
  z <- stats::qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE)
  n <- nrow(dataset)
  if (full_denom) {
    n_Ic <- n
    n_It <- n
    n_IcMc <- n
    n_IcMt <- n
    n_ItMc <- n
    n_ItMt <- n
  } else {
    n_Ic <- nrow(subset(dataset, induction == "Ic"))
    n_It <- nrow(subset(dataset, induction == "It"))
    n_IcMc <- nrow(subset(dataset, induction == "Ic" & maintenance == "Mc"))
    n_IcMt <- nrow(subset(dataset, induction == "Ic" & maintenance == "Mt"))
    n_ItMc <- nrow(subset(dataset, induction == "It" & maintenance == "Mc"))
    n_ItMt <- nrow(subset(dataset, induction == "It" & maintenance == "Mt"))
  }
  theta_Ic_ci <- c(
    paramters$theta_Ic - (z / sqrt(fisher_info$theta_Ic_info * n_Ic)),
    paramters$theta_Ic + (z / sqrt(fisher_info$theta_Ic_info * n_Ic))
  )
  theta_It_ci <- c(
    paramters$theta_It - (z / sqrt(fisher_info$theta_It_info * n_It)),
    paramters$theta_It + (z / sqrt(fisher_info$theta_It_info * n_It))
  )
  lambda_Ic_nr_ci <- c(
    paramters$lambda_Ic_nr - (z / sqrt(fisher_info$lambda_Ic_nr_info * n_Ic)),
    paramters$lambda_Ic_nr + (z / sqrt(fisher_info$lambda_Ic_nr_info * n_Ic))
  )
  lambda_It_nr_ci <- c(
    paramters$lambda_It_nr - (z / sqrt(fisher_info$lambda_It_nr_info * n_It)),
    paramters$lambda_It_nr + (z / sqrt(fisher_info$lambda_It_nr_info * n_It))
  )
  lambda_IcMc_r_ci <- c(
    paramters$lambda_IcMc_r - (z / sqrt(fisher_info$lambda_IcMc_r_info * n_IcMc)),
    paramters$lambda_IcMc_r + (z / sqrt(fisher_info$lambda_IcMc_r_info * n_IcMc))
  )
  lambda_IcMc_nr_ci <- c(
    paramters$lambda_IcMc_nr - (z / sqrt(fisher_info$lambda_IcMc_nr_info * n_IcMc)),
    paramters$lambda_IcMc_nr + (z / sqrt(fisher_info$lambda_IcMc_nr_info * n_IcMc))
  )
  lambda_IcMt_r_ci <- c(
    paramters$lambda_IcMt_r - (z / sqrt(fisher_info$lambda_IcMt_r_info * n_IcMt)),
    paramters$lambda_IcMt_r + (z / sqrt(fisher_info$lambda_IcMt_r_info * n_IcMt))
  )
  lambda_IcMt_nr_ci <- c(
    paramters$lambda_IcMt_nr - (z / sqrt(fisher_info$lambda_IcMt_nr_info * n_IcMt)),
    paramters$lambda_IcMt_nr + (z / sqrt(fisher_info$lambda_IcMt_nr_info * n_IcMt))
  )
  lambda_ItMc_r_ci <- c(
    paramters$lambda_ItMc_r - (z / sqrt(fisher_info$lambda_ItMc_r_info * n_ItMc)),
    paramters$lambda_ItMc_r + (z / sqrt(fisher_info$lambda_ItMc_r_info * n_ItMc))
  )
  lambda_ItMc_nr_ci <- c(
    paramters$lambda_ItMc_nr - (z / sqrt(fisher_info$lambda_ItMc_nr_info * n_ItMc)),
    paramters$lambda_ItMc_nr + (z / sqrt(fisher_info$lambda_ItMc_nr_info * n_ItMc))
  )
  lambda_ItMt_r_ci <- c(
    paramters$lambda_ItMt_r - (z / sqrt(fisher_info$lambda_ItMt_r_info * n_ItMt)),
    paramters$lambda_ItMt_r + (z / sqrt(fisher_info$lambda_ItMt_r_info * n_ItMt))
  )
  lambda_ItMt_nr_ci <- c(
    paramters$lambda_ItMt_nr - (z / sqrt(fisher_info$lambda_ItMt_nr_info * n_ItMt)),
    paramters$lambda_ItMt_nr + (z / sqrt(fisher_info$lambda_ItMt_nr_info * n_ItMt))
  )
  result <- list(
    theta_Ic_ci = theta_Ic_ci,
    theta_It_ci = theta_It_ci,
    lambda_Ic_nr_ci = lambda_Ic_nr_ci,
    lambda_It_nr_ci = lambda_It_nr_ci,
    lambda_IcMc_r_ci = lambda_IcMc_r_ci,
    lambda_IcMc_nr_ci = lambda_IcMc_nr_ci,
    lambda_IcMt_r_ci = lambda_IcMt_r_ci,
    lambda_IcMt_nr_ci = lambda_IcMt_nr_ci,
    lambda_ItMc_r_ci = lambda_ItMc_r_ci,
    lambda_ItMc_nr_ci = lambda_ItMc_nr_ci,
    lambda_ItMt_r_ci = lambda_ItMt_r_ci,
    lambda_ItMt_nr_ci = lambda_ItMt_nr_ci
  )
  return (result)
}

#' Parameter estimation using EM algorithm
#'
#' This function is to estimate parameters using EM algorithm.
#' The return value is list.
#' One dataset from generate_scenario is used in estimation.
#'
#' @param dataset The dataset of one scenario from function generate_scenario.
#'
#' @param theta_Ic_init An initial value of proportion of responders in induction control.
#'
#' @param theta_It_init An initial value of proportion of responders in induction treatment.
#'
#' @param lambda_Ic_nr_init An initial value of hazard of non-responders in induction control.
#'
#' @param lambda_It_nr_init An initial value of hazard of non-responders in induction treatment.
#'
#' @param lambda_IcMc_r_init An initial value of hazard of responders in maintenance control when induction control.
#'
#' @param lambda_IcMc_nr_init An initial value of hazard of non-responders in maintenance control when induction control.
#'
#' @param lambda_ItMc_r_init An initial value of hazard of responders in maintenance treatment when induction control.
#'
#' @param lambda_ItMc_nr_init An initial value of hazard of non-responders in maintenance treatment when induction control.
#'
#' @param lambda_IcMt_r_init An initial value of hazard of responders in maintenance control when induction treatment.
#'
#' @param lambda_IcMt_nr_init An initial value of hazard of non-responders in maintenance control when induction treatment.
#'
#' @param lambda_ItMt_r_init An initial value of hazard of responders in maintenance treatment when induction treatment.
#'
#' @param lambda_ItMt_nr_init An initial value of hazard of non-responders in maintenance treatment when induction treatment.
#'
#' @param max_iter The maximum number of iteration in each EM algorithm.
#'
#' @param eps The threshold of difference log-likelihoods to stop EM algorithm.
#'
#' @param option If NULL, implement default EM algorithm.
#' If "fix_theta", implement EM algorithm being theta fixed.
#'
#' @export
estimateEM <- function (
  dataset,
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
  # fixed data
  X_IcMc <- subset(
    dataset, induction == "Ic" & maintenance == "Mc",
    select = c(time, cens, t_judge, is_induction)
  )
  X_ItMc <- subset(
    dataset, induction == "It" & maintenance == "Mc",
    select = c(time, cens, t_judge, is_induction)
  )
  X_IcMt <- subset(
    dataset, induction == "Ic" & maintenance == "Mt",
    select = c(time, cens, t_judge, is_induction)
  )
  X_ItMt <- subset(
    dataset, induction == "It" & maintenance == "Mt",
    select = c(time, cens, t_judge, is_induction)
  )
  true_theta_Ic <- subset(
    dataset, induction == "Ic", select = theta
  )[1,1]
  true_theta_It <- subset(
    dataset, induction == "It", select = theta
  )[1,1]

  # parameters to be updated
  theta_Ic <- theta_Ic_init
  theta_It <- theta_It_init
  lambda_Ic_nr <- lambda_Ic_nr_init
  lambda_It_nr <- lambda_It_nr_init
  lambda_IcMc_r <- lambda_IcMc_r_init
  lambda_IcMc_nr <- lambda_IcMc_nr_init
  lambda_ItMc_r <- lambda_ItMc_r_init
  lambda_ItMc_nr <- lambda_ItMc_nr_init
  lambda_IcMt_r <- lambda_IcMt_r_init
  lambda_IcMt_nr <- lambda_IcMt_nr_init
  lambda_ItMt_r <- lambda_ItMt_r_init
  lambda_ItMt_nr <- lambda_ItMt_nr_init

  pi_IcMc <- NULL
  pi_IcMt <- NULL
  pi_ItMc <- NULL
  pi_ItMt <- NULL


  # define functions to update parameters and to calculate log-likelihood
  if (is.null(option)) {
    # default EM algorithm
    Estep <- function () {
      pi_IcMc <<- calc_responsibility(
        X_IM = X_IcMc,
        theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
      )
      pi_IcMt <<- calc_responsibility(
        X_IM = X_IcMt,
        theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
      )
      pi_ItMc <<- calc_responsibility(
        X_IM = X_ItMc,
        theta = theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
      )
      pi_ItMt <<- calc_responsibility(
        X_IM = X_ItMt,
        theta = theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
      )
    }

    Mstep <- function () {
      theta_Ic <<- update_theta(
        X_IMc = X_IcMc, X_IMt = X_IcMt,
        pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
      )
      theta_It <<- update_theta(
        X_IMc = X_ItMc, X_IMt = X_ItMt,
        pi_IMc = pi_ItMc, pi_IMt = pi_ItMt
      )

      lambda_Ic_nr <<- update_lambda_I(
        X_IMc = X_IcMc, X_IMt = X_IcMt,
        pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
      )
      lambda_It_nr <<- update_lambda_I(
        X_IMc = X_ItMc, X_IMt = X_ItMt,
        pi_IMc = pi_ItMc, pi_IMt = pi_ItMt
      )

      lambda_IcMc <- update_lambda_IM(X_IM = X_IcMc, pi_IM = pi_IcMc)
      lambda_ItMc <- update_lambda_IM(X_IM = X_ItMc, pi_IM = pi_ItMc)
      lambda_IcMt <- update_lambda_IM(X_IM = X_IcMt, pi_IM = pi_IcMt)
      lambda_ItMt <- update_lambda_IM(X_IM = X_ItMt, pi_IM = pi_ItMt)
      lambda_IcMc_r <<- lambda_IcMc$r
      lambda_IcMc_nr <<- lambda_IcMc$nr
      lambda_ItMc_r <<- lambda_ItMc$r
      lambda_ItMc_nr <<- lambda_ItMc$nr
      lambda_IcMt_r <<- lambda_IcMt$r
      lambda_IcMt_nr <<- lambda_IcMt$nr
      lambda_ItMt_r <<- lambda_ItMt$r
      lambda_ItMt_nr <<- lambda_ItMt$nr
    }
  } else if (option == "fix_theta") {
    # method to fix theta
    theta_Ic <- true_theta_Ic
    theta_It <- true_theta_It

    Estep <- function () {
      pi_IcMc <<- calc_responsibility(
        X_IM = X_IcMc,
        theta = true_theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
      )
      pi_IcMt <<- calc_responsibility(
        X_IM = X_IcMt,
        theta = true_theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
      )
      pi_ItMc <<- calc_responsibility(
        X_IM = X_ItMc,
        theta = true_theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
      )
      pi_ItMt <<- calc_responsibility(
        X_IM = X_ItMt,
        theta = true_theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
      )
    }

    Mstep <- function () {
      lambda_Ic_nr <<- update_lambda_I(
        X_IMc = X_IcMc, X_IMt = X_IcMt,
        pi_IMc = pi_IcMc, pi_IMt = pi_IcMt
      )
      lambda_It_nr <<- update_lambda_I(
        X_IMc = X_ItMc, X_IMt = X_ItMt,
        pi_IMc = pi_ItMc, pi_IMt = pi_ItMt
      )

      lambda_IcMc <- update_lambda_IM(X_IM = X_IcMc, pi_IM = pi_IcMc)
      lambda_ItMc <- update_lambda_IM(X_IM = X_ItMc, pi_IM = pi_ItMc)
      lambda_IcMt <- update_lambda_IM(X_IM = X_IcMt, pi_IM = pi_IcMt)
      lambda_ItMt <- update_lambda_IM(X_IM = X_ItMt, pi_IM = pi_ItMt)
      lambda_IcMc_r <<- lambda_IcMc$r
      lambda_IcMc_nr <<- lambda_IcMc$nr
      lambda_ItMc_r <<- lambda_ItMc$r
      lambda_ItMc_nr <<- lambda_ItMc$nr
      lambda_IcMt_r <<- lambda_IcMt$r
      lambda_IcMt_nr <<- lambda_IcMt$nr
      lambda_ItMt_r <<- lambda_ItMt$r
      lambda_ItMt_nr <<- lambda_ItMt$nr
    }
  } else {
    stop("option parameter is invalid.")
  }

  calc_loglikelihood <- function () {
    loglikelihood <- (
      calc_loglikelihood_stage1(
        X_IMc = X_IcMc, X_IMt = X_IcMt,
        theta = theta_Ic, lambda_I_nr = lambda_Ic_nr
      ) +
      calc_loglikelihood_stage1(
        X_IMc = X_ItMc, X_IMt = X_ItMt,
        theta = theta_It, lambda_I_nr = lambda_It_nr
      ) +
      calc_loglikelihood_stage2(
        X_IM = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
      ) +
      calc_loglikelihood_stage2(
        X_IM = X_ItMc, theta = theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
      ) +
      calc_loglikelihood_stage2(
        X_IM = X_IcMt, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
        lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
      ) +
      calc_loglikelihood_stage2(
        X_IM = X_ItMt, theta = theta_It, lambda_I_nr = lambda_It_nr,
        lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
      )
    )
    return (loglikelihood)
  }

  loglikelihoods <- rep(NULL, max_iter)
  last_parameters <- list()

  # implement EM algorithm
  for (i in 1:max_iter) {
    Estep()
    Mstep()
    loglikelihoods[i] <- calc_loglikelihood()
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
      lambda_IcMc_r = lambda_IcMc_r,
      lambda_IcMc_nr = lambda_IcMc_nr,
      lambda_ItMc_r = lambda_ItMc_r,
      lambda_ItMc_nr = lambda_ItMc_nr,
      lambda_IcMt_r = lambda_IcMt_r,
      lambda_IcMt_nr = lambda_IcMt_nr,
      lambda_ItMt_r = lambda_ItMt_r,
      lambda_ItMt_nr = lambda_ItMt_nr
    )
    if (i > 1) {
      if ((loglikelihoods[i] - loglikelihoods[i-1]) < 0) {
        loglikelihoods <- loglikelihoods[1:i-1]
        parameters <- last_parameters
        break
      }
      if (abs(loglikelihoods[i] - loglikelihoods[i-1]) <= eps) { break }
    }
    last_parameters <- parameters
  }

  loglikelihoods <- loglikelihoods[!is.null(loglikelihoods)]
  step <- length(loglikelihoods) - 1
  fisher_information <- calc_Fisher_information(dataset, parameters)

  return (c(
    list(
      step = step,
      loglikelihood = loglikelihoods
    ),
    parameters,
    calc_Fisher_information(dataset, parameters)
  ))
}

#' Parameter estimation in the form of data.frame
#'
#' This function is to estimate parameters using EM algorithm,
#' and is the wrapper of estimateEM function.
#' The return value is data.frame.
#' One dataset from generate_scenario is used in estimation.
#'
#' @param dataset The dataset of one scenario from function generate_scenario.
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
estimateEM_as_frame <- function (
    dataset,
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
  n_IcMc <- nrow(subset(dataset, induction == "Ic" & maintenance == "Mc"))
  n_ItMc <- nrow(subset(dataset, induction == "It" & maintenance == "Mc"))
  n_IcMt <- nrow(subset(dataset, induction == "Ic" & maintenance == "Mt"))
  n_ItMt <- nrow(subset(dataset, induction == "It" & maintenance == "Mt"))
  t_judge <- dataset$t_judge[1]

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
    lambda_ItMt_nr_info = e$lambda_ItMt_nr_info,
    # loglikelihood at the last step
    last_loglikelihood = e$loglikelihood[e$step + 1],
    # sample size
    n_IcMc = n_IcMc,
    n_ItMc = n_ItMc,
    n_IcMt = n_IcMt,
    n_ItMt = n_ItMt,
    t_judge = t_judge
  )
  return (estimator)
}

#' Parameter estimation by sequential execution
#'
#' This function is to estimate parameters by sequential execution.
#' Full dataset from generate_scenario is used in estimation.
#'
#' @param dataset The datasets of scenarios from function generate_scenario.
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
estimate_sequentially <- function (
    dataset,
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
  estimate_func <- function (data) {
    target_id <- data$id[1]
    if (verbose) {
      cat(target_id, fill = TRUE)
    }
    result_frame <- estimateEM_as_frame(
      data,
      theta_Ic_init = theta_Ic_init[target_id],
      theta_It_init = theta_Ic_init[target_id],
      lambda_Ic_nr_init = lambda_Ic_nr_init[target_id],
      lambda_It_nr_init = lambda_It_nr_init[target_id],
      lambda_IcMc_r_init = lambda_IcMc_r_init[target_id],
      lambda_IcMc_nr_init = lambda_IcMc_nr_init[target_id],
      lambda_ItMc_r_init = lambda_ItMc_r_init[target_id],
      lambda_ItMc_nr_init = lambda_ItMc_nr_init[target_id],
      lambda_IcMt_r_init = lambda_IcMt_r_init[target_id],
      lambda_IcMt_nr_init = lambda_IcMt_nr_init[target_id],
      lambda_ItMt_r_init = lambda_ItMt_r_init[target_id],
      lambda_ItMt_nr_init = lambda_ItMt_nr_init[target_id],
      max_iter = max_iter,
      eps = eps,
      option = option
    )
    cbind(list(id = target_id), result_frame)
  }

  # execution
  estimators <- lapply(
    split_dataset, estimate_func
  )
  result <- Reduce(rbind, estimators)
  result
}

replicate_dataset <- function (
  original_dataset,
  replicate_num = 50L
) {
  renumber_id_func <- function (new_id) {
    rep_dataset <- original_dataset
    rep_dataset$id <- new_id
    rep_dataset
  }
  replicate_ids <- 1:replicate_num
  replicated_list <- lapply(replicate_ids, renumber_id_func)
  Reduce(rbind, replicated_list)
}

#' Parameter estimation by sequential execution to get max loglikelihood
#'
#' This function is to estimate parameters by sequential execution.
#' Full dataset from generate_scenario is used in estimation.
#'
#' @param dataset The datasets of scenarios from function generate_scenario.
#'
#' @param verbose If TRUE, show progress.
#'
#' @param maximize_num The number of iteration.
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
estimate_sequentially_max_likelihood <- function (
  dataset,
  verbose = TRUE,
  maximize_num = 100L,
  theta_Ic_init = rep(0.65, 100L),
  theta_It_init = rep(0.75, 100L),
  lambda_Ic_nr_init = rep(0.10, 100L),
  lambda_It_nr_init = rep(0.05, 100L),
  lambda_IcMc_r_init = rep(0.10, 100L),
  lambda_IcMc_nr_init = rep(0.20, 100L),
  lambda_ItMc_r_init = rep(0.08, 100L),
  lambda_ItMc_nr_init = rep(0.16, 100L),
  lambda_IcMt_r_init = rep(0.08, 100L),
  lambda_IcMt_nr_init = rep(0.06, 100L),
  lambda_ItMt_r_init = rep(0.04, 100L),
  lambda_ItMt_nr_init = rep(0.08, 100L),
  max_iter = 5000L,
  eps = 1e-5,
  option = NULL
) {
  target_ids <- unique(dataset$id)
  split_dataset <- lapply(
    target_ids,
    function (target_id) {dataset[dataset$id == target_id,]}
  )
  estimate_func <- function (dataset) {
    target_id <- dataset$id[1]
    if (verbose) {
      cat(target_id, fill = TRUE)
    }

    # achieve result of the max loglikelihood
    replicated_datasets <- replicate_dataset(dataset, replicate_num = maximize_num)
    repeat_estimate <- estimate_sequentially(
      replicated_datasets,
      verbose = verbose,
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
    max_loglikelihood_index <- which.max(repeat_estimate$last_loglikelihood)
    max_loglikelihood_row <- subset(repeat_estimate[max_loglikelihood_index,], select = -id)

    cbind(list(id = target_id), max_loglikelihood_row)
  }

  # execution
  estimators <- lapply(
    split_dataset, estimate_func
  )
  result <- Reduce(rbind, estimators)
  result
}
