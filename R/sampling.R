rexp_range <- function (
  n,
  lambda,
  range = c(0, 1)
) {
  u <- runif(n, min = range[1], max = range[2])
  t <- -log(u)/lambda
  return (t)
}

sample_mixture_exp <- function (
  n,
  theta,
  lambda1,
  lambda2,
  range1 = c(0, 1),
  range2 = c(0, 1)
) {
  if (lambda2 <= lambda1) {
    stop("lambda2 should be greater than lambda1.")
  }

  u <- runif(n)
  is_less <- ifelse(u < theta, 1L, 0L)
  n1 <- sum(is_less)
  n2 <- n - n1
  if (lambda1 == 0) {
    t <- rexp_range(n, lambda2, range2)
  } else {
    t <- c(
      rexp_range(n1, lambda1, range1),
      rexp_range(n2, lambda2, range2)
    )
  }

  return (t)
}

#' @export
sample_tsrd_dataset <- function (
  n,
  t_judge,
  theta,
  lambda_I_r,
  lambda_I_nr,
  lambda_IM_r,
  lambda_IM_nr
) {
  S_r_judge <- exp(-lambda_I_r * t_judge)
  S_nr_judge <- exp(-lambda_I_nr * t_judge)
  S_judge <- theta * S_r_judge + (1 - theta) * S_nr_judge
  S <- runif(n)
  is_induction <- ifelse(S >= S_judge, 1L, 0L)
  n_I <- sum(is_induction)
  n_IM <- n - n_I
  t_I <- sample_mixture_exp(
    n_I, theta, lambda_I_r, lambda_I_nr,
    range1 = c(S_r_judge, 1), range2 = c(S_nr_judge, 1)
  )
  t_IM <- sample_mixture_exp(
    n_IM, theta, lambda_IM_r, lambda_IM_nr
  )
  t <- c(t_I, (t_IM + t_judge))
  return (t)
}
