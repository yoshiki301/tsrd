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

#' Parameter estimation using EM algorithm
#'
#' This function is to estimate parameters using EM algorithm.
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

  return (c(
    list(
      step = step,
      loglikelihood = loglikelihoods
    ),
    parameters
  ))
}
