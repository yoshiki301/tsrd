f <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  f_r <- a * 0.0 + (1 - a) * stats::dexp(t - t_judge, lambda_IM_r)
  f_nr <- a * stats::dexp(t, lambda_I_nr) + (1 - a) * c_nr * stats::dexp(t - t_judge, lambda_IM_nr)

  theta * f_r + (1.0 - theta) * f_nr
}

f_p_theta <- function (
  t,
  t_judge,
  a,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  f_r <- a * 0.0 + (1 - a) * stats::dexp(t - t_judge, lambda_IM_r)
  f_nr <- a * stats::dexp(t, lambda_I_nr) + (1 - a) * c_nr * stats::dexp(t - t_judge, lambda_IM_nr)

  f_r - f_nr
}

f_p_lambda_I_nr <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_I_nr,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  induction <- (1 - lambda_I_nr * t) * (1.0 - stats::pexp(t, lambda_I_nr))
  maintenance <- t_judge * c_nr * stats::dexp(t - t_judge, lambda_IM_nr)

  -(1.0 - theta) * (a * induction + (1 - a) * maintenance)
}

f_p_lambda_IM_r <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_IM_r
) {
  theta * (1 - a) * (1 - lambda_IM_r * (t - t_judge)) * (1.0 - stats::pexp(t - t_judge, lambda_IM_r))
}

f_p_lambda_IM_nr <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_I_nr,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  (1.0 - theta) * (1 - a) * (1 - lambda_IM_nr * (t - t_judge)) * c_nr * (1.0 - stats::pexp(t - t_judge, lambda_IM_nr))
}

S <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  S_r <- a * 1.0 + (1 - a) * (1.0 - stats::pexp(t - t_judge, lambda_IM_r))
  S_nr <- a * (1.0 - stats::pexp(t, lambda_I_nr)) + (1 - a) * c_nr * (1.0 - stats::pexp(t - t_judge, lambda_IM_nr))

  theta * S_r + (1.0 - theta) * S_nr
}

S_p_theta <- function (
  t,
  t_judge,
  a,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  S_r <- a * 1.0 + (1 - a) * (1.0 - stats::pexp(t - t_judge, lambda_IM_r))
  S_nr <- a * (1.0 - stats::pexp(t, lambda_I_nr)) + (1 - a) * c_nr * (1.0 - stats::pexp(t - t_judge, lambda_IM_nr))

  S_r - S_nr
}

S_p_lambda_I_nr <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_I_nr,
  lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  induction <- t * (1.0 - stats::pexp(t, lambda_I_nr))
  maintenance <- t_judge * c_nr * (1.0 - stats::pexp(t - t_judge, lambda_IM_nr))

  -(1.0 - theta) * (a * induction + (1 - a) * maintenance)
}

S_p_lambda_IM_r <- function (
  t,
  t_judge,
  a,
  theta,
  lambda_IM_r
) {
  theta * (1 - a) * (t - t_judge) * (1.0 - stats::pexp(t - t_judge, lambda_IM_r))
}

S_p_lambda_IM_nr <- function (
    t,
    t_judge,
    a,
    theta,
    lambda_I_nr,
    lambda_IM_nr
) {
  c_nr <- 1.0 - stats::pexp(t_judge, lambda_I_nr)

  (1.0 - theta) * (1 - a) * (t - t_judge) * c_nr * (1.0 - stats::pexp(t - t_judge, lambda_IM_nr))
}

lp_theta <- function (
  X,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  X$cens * (
    f_p_theta(
      t = X$time, t_judge = X$t_judge, a = X$is_induction,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    ) / f(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  ) + (1 - X$cens) * (
    S_p_theta(
      t = X$time, t_judge = X$t_judge, a = X$is_induction,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    ) / S(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  )
}

lp_lambda_I_nr <- function (
  X,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  X$cens * (
    f_p_lambda_I_nr(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_nr = lambda_IM_nr
    ) / f(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  ) + (1 - X$cens) * (
    S_p_lambda_I_nr(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_nr = lambda_IM_nr
    ) / S(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  )
}

lp_lambda_IM_r <- function (
  X,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  X$cens * (
    f_p_lambda_IM_r(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_IM_r = lambda_IM_r
    ) / f(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  ) + (1 - X$cens) * (
    S_p_lambda_IM_r(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_IM_r = lambda_IM_r
    ) / S(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  )
}

lp_lambda_IM_nr <- function (
  X,
  theta,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  X$cens * (
    f_p_lambda_IM_nr(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_nr = lambda_IM_nr
    ) / f(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  ) + (1 - X$cens) * (
    S_p_lambda_IM_nr(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_nr = lambda_IM_nr
    ) / S(
      t = X$time, t_judge = X$t_judge, a = X$is_induction, theta = theta,
      lambda_I_nr = lambda_I_nr, lambda_IM_r = lambda_IM_r, lambda_IM_nr = lambda_IM_nr
    )
  )
}
