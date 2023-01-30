adjust_censored_time <- function (
  time,
  cens,
  seed = 42L
) {
  set.seed(seed = seed)
  if (length(time) != length(cens)) {
    stop("The lengths of time and cens should be equal.")
  }

  censored_time <- time[cens == 0]
  time[cens == 0] <- stats::runif(length(censored_time), 0.0, censored_time)
  return (time)
}

#' Generate datasets in one scenario
#'
#' This function is to generate datasets with true parameters in each treatment group.
#' Each dataset consists of survival time, censoring status, and treatment assignment.
#'
#' @param sim_num The number of simulations.
#'
#' @param t_judge The time to judge whether to be responder or not.
#'
#' @param t_max The time of the end of trial.
#'
#' @param censor_rate The proportion of censoring.
#'
#' @param IcMc_size The sample size of induction control and maintenance control.
#'
#' @param ItMc_size The sample size of induction treatment and maintenance control.
#'
#' @param IcMt_size The sample size of induction control and maintenance treatment.
#'
#' @param ItMt_size The sample size of induction treatment and maintenance treatment.
#'
#' @param theta_Ic The proportion of responders in induction control.
#'
#' @param theta_It The proportion of responders in induction treatment.
#'
#' @param lambda_Ic_r The hazard of responders in induction control. This should be set to 0.
#'
#' @param lambda_Ic_nr The hazard of non-responders in induction control.
#'
#' @param lambda_It_r The hazard of responders in induction treatment. This should be set to 0.
#'
#' @param lambda_It_nr The hazard of non-responders in induction treatment.
#'
#' @param lambda_IcMc_r The hazard of responders in maintenance control when induction control.
#'
#' @param lambda_IcMc_nr The hazard of non-responders in maintenance control when induction control.
#'
#' @param lambda_ItMc_r The hazard of responders in maintenance treatment when induction control.
#'
#' @param lambda_ItMc_nr The hazard of non-responders in maintenance treatment when induction control.
#'
#' @param lambda_IcMt_r The hazard of responders in maintenance control when induction treatment.
#'
#' @param lambda_IcMt_nr The hazard of non-responders in maintenance control when induction treatment.
#'
#' @param lambda_ItMt_r The hazard of responders in maintenance treatment when induction treatment.
#'
#' @param lambda_ItMt_nr The hazard of non-responders in maintenance treatment when induction treatment.
#'
#' @param seed The value of random seed.
#'
#' @export
generate_scenario <- function (
  sim_num = 1000L,
  t_judge = 6,
  t_max = 60,
  censor_rate = 0.0,
  IcMc_size = 1000L,
  ItMc_size = 1000L,
  IcMt_size = 1000L,
  ItMt_size = 1000L,
  theta_Ic = 0.65,
  theta_It = 0.75,
  lambda_Ic_r = 0.0,
  lambda_Ic_nr = 0.10,
  lambda_It_r = 0.0,
  lambda_It_nr = 0.05,
  lambda_IcMc_r = 0.10,
  lambda_IcMc_nr = 0.20,
  lambda_ItMc_r = 0.08,
  lambda_ItMc_nr = 0.16,
  lambda_IcMt_r = 0.08,
  lambda_IcMt_nr = 0.16,
  lambda_ItMt_r = 0.04,
  lambda_ItMt_nr = 0.08,
  seed = 42L
) {
  seeds <- (1L:sim_num) + seed

  generate_func <- function (sim_id) {
    set.seed(seeds[sim_id])

    # sample random variable of time
    # This is ideal time, observation is censored or not
    t_IcMc <- sample_tsrd_dataset(
      IcMc_size, t_judge, theta_Ic,
      lambda_Ic_r, lambda_Ic_nr,
      lambda_IcMc_r, lambda_IcMc_nr
    )
    t_ItMc <- sample_tsrd_dataset(
      ItMc_size, t_judge, theta_It,
      lambda_It_r, lambda_It_nr,
      lambda_ItMc_r, lambda_ItMc_nr
    )
    t_IcMt <- sample_tsrd_dataset(
      IcMt_size, t_judge, theta_Ic,
      lambda_Ic_r, lambda_Ic_nr,
      lambda_IcMt_r, lambda_IcMt_nr
    )
    t_ItMt <- sample_tsrd_dataset(
      ItMt_size, t_judge, theta_It,
      lambda_It_r, lambda_It_nr,
      lambda_ItMt_r, lambda_ItMt_nr
    )

    # set status of censored or not
    cens_IcMc <- stats::rbinom(IcMc_size, 1L, 1.0 - censor_rate)
    cens_ItMc <- stats::rbinom(ItMc_size, 1L, 1.0 - censor_rate)
    cens_IcMt <- stats::rbinom(IcMt_size, 1L, 1.0 - censor_rate)
    cens_ItMt <- stats::rbinom(ItMt_size, 1L, 1.0 - censor_rate)

    # calculate observed time (event or censoring)
    obs_IcMc <- adjust_censored_time(t_IcMc, cens_IcMc)
    obs_ItMc <- adjust_censored_time(t_ItMc, cens_ItMc)
    obs_IcMt <- adjust_censored_time(t_IcMt, cens_IcMt)
    obs_ItMt <- adjust_censored_time(t_ItMt, cens_ItMt)

    # set censoring at the end of trial
    obs_IcMc <- ifelse(obs_IcMc <= t_max, obs_IcMc, t_max)
    cens_IcMc <- ifelse(obs_IcMc < t_max, cens_IcMc, 0L)
    obs_ItMc <- ifelse(obs_ItMc <= t_max, obs_ItMc, t_max)
    cens_ItMc <- ifelse(obs_ItMc < t_max, cens_ItMc, 0L)
    obs_IcMt <- ifelse(obs_IcMt <= t_max, obs_IcMt, t_max)
    cens_IcMt <- ifelse(obs_IcMt < t_max, cens_IcMt, 0L)
    obs_ItMt <- ifelse(obs_ItMt <= t_max, obs_ItMt, t_max)
    cens_ItMt <- ifelse(obs_ItMt < t_max, cens_ItMt, 0L)

    # make a scenario data as data.frame
    data_IcMc <- data.frame(
      id = sim_id,
      time = obs_IcMc,
      cens = cens_IcMc,
      t_judge = t_judge,
      t_max = t_max,
      is_induction = ifelse(obs_IcMc <= t_judge, 1L, 0L),
      induction = "Ic",
      maintenance = "Mc",
      theta = theta_Ic,
      lambda_I_r = lambda_Ic_r,
      lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMc_r,
      lambda_IM_nr = lambda_IcMc_nr
    )
    data_ItMc <- data.frame(
      id = sim_id,
      time = obs_ItMc,
      cens = cens_ItMc,
      t_judge = t_judge,
      t_max = t_max,
      is_induction = ifelse(obs_ItMc <= t_judge, 1L, 0L),
      induction = "It",
      maintenance = "Mc",
      theta = theta_It,
      lambda_I_r = lambda_It_r,
      lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMc_r,
      lambda_IM_nr = lambda_ItMc_nr
    )
    data_IcMt <- data.frame(
      id = sim_id,
      time = obs_IcMt,
      cens = cens_IcMt,
      t_judge = t_judge,
      t_max = t_max,
      is_induction = ifelse(obs_IcMt <= t_judge, 1L, 0L),
      induction = "Ic",
      maintenance = "Mt",
      theta = theta_Ic,
      lambda_I_r = lambda_Ic_r,
      lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMt_r,
      lambda_IM_nr = lambda_IcMt_nr
    )
    data_ItMt <- data.frame(
      id = sim_id,
      time = obs_ItMt,
      cens = cens_ItMt,
      t_judge = t_judge,
      t_max = t_max,
      is_induction = ifelse(obs_ItMt <= t_judge, 1L, 0L),
      induction = "It",
      maintenance = "Mt",
      theta = theta_It,
      lambda_I_r = lambda_It_r,
      lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMt_r,
      lambda_IM_nr = lambda_ItMt_nr
    )

    data <- rbind(data_IcMc, data_ItMc, data_IcMt, data_ItMt)
    return (data)
  }

  datasets <- lapply(1L:sim_num, generate_func)
  scenario <- Reduce(rbind, datasets)
  return (scenario)
}
