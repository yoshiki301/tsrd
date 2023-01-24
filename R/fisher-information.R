loglikelihood_grad_vector <- function (
  dataset,
  estimators
) {
  # fixed data
  X_IcMc <- subset(
    dataset, induction == "Ic" & maintenance == "Mc",
    select = c(time, t_judge, cens, is_induction)
  )
  X_ItMc <- subset(
    dataset, induction == "It" & maintenance == "Mc",
    select = c(time, t_judge, cens, is_induction)
  )
  X_IcMt <- subset(
    dataset, induction == "Ic" & maintenance == "Mt",
    select = c(time, t_judge, cens, is_induction)
  )
  X_ItMt <- subset(
    dataset, induction == "It" & maintenance == "Mt",
    select = c(time, t_judge, cens, is_induction)
  )

  # estimators of parameter
  theta_Ic <- estimators$theta_Ic
  theta_It <- estimators$theta_It
  lambda_Ic_nr <- estimators$lambda_Ic_nr
  lambda_It_nr <- estimators$lambda_It_nr
  lambda_IcMc_r <- estimators$lambda_IcMc_r
  lambda_IcMc_nr <- estimators$lambda_IcMc_nr
  lambda_ItMc_r <- estimators$lambda_ItMc_r
  lambda_ItMc_nr <- estimators$lambda_ItMc_nr
  lambda_IcMt_r <- estimators$lambda_IcMt_r
  lambda_IcMt_nr <- estimators$lambda_IcMt_nr
  lambda_ItMt_r <- estimators$lambda_ItMt_r
  lambda_ItMt_nr <- estimators$lambda_ItMt_nr

  # gradient of loglikelihood
  lp_theta_Ic <- sum(
    lp_theta(
      X = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
    ),
    lp_theta(
      X = X_IcMt, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
    )
  )
  lp_theta_It <- sum(
    lp_theta(
      X = X_ItMc, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
    ),
    lp_theta(
      X = X_ItMt, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
    )
  )
  lp_lambda_Ic_nr <- sum(
    lp_lambda_I_nr(
      X = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
    ),
    lp_lambda_I_nr(
      X = X_IcMt, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
    )
  )
  lp_lambda_It_nr <- sum(
    lp_lambda_I_nr(
      X = X_ItMc, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
    ),
    lp_lambda_I_nr(
      X = X_ItMt, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
    )
  )
  lp_lambda_IcMc_r <- sum(
    lp_lambda_IM_r(
      X = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
    )
  )
  lp_lambda_IcMc_nr <- sum(
    lp_lambda_IM_nr(
      X = X_IcMc, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMc_r, lambda_IM_nr = lambda_IcMc_nr
    )
  )
  lp_lambda_ItMc_r <- sum(
    lp_lambda_IM_r(
      X = X_ItMc, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
    )
  )
  lp_lambda_ItMc_nr <- sum(
    lp_lambda_IM_nr(
      X = X_ItMc, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMc_r, lambda_IM_nr = lambda_ItMc_nr
    )
  )
  lp_lambda_IcMt_r <- sum(
    lp_lambda_IM_r(
      X = X_IcMt, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
    )
  )
  lp_lambda_IcMt_nr <- sum(
    lp_lambda_IM_nr(
      X = X_IcMt, theta = theta_Ic, lambda_I_nr = lambda_Ic_nr,
      lambda_IM_r = lambda_IcMt_r, lambda_IM_nr = lambda_IcMt_nr
    )
  )
  lp_lambda_ItMt_r <- sum(
    lp_lambda_IM_r(
      X = X_ItMt, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
    )
  )
  lp_lambda_ItMt_nr <- sum(
    lp_lambda_IM_nr(
      X = X_ItMt, theta = theta_It, lambda_I_nr = lambda_It_nr,
      lambda_IM_r = lambda_ItMt_r, lambda_IM_nr = lambda_ItMt_nr
    )
  )

  c(
    lp_theta_Ic,
    lp_theta_It,
    lp_lambda_Ic_nr,
    lp_lambda_It_nr,
    lp_lambda_IcMc_r,
    lp_lambda_IcMc_nr,
    lp_lambda_ItMc_r,
    lp_lambda_ItMc_nr,
    lp_lambda_IcMt_r,
    lp_lambda_IcMt_nr,
    lp_lambda_ItMt_r,
    lp_lambda_ItMt_nr
  )
}

generate_pseudo_dataset <- function (
  original_dataset,
  estimators,
  num = 100L,
  censor_rate = 0.0,
  seed = 42L
) {
  n_IcMc <- nrow(subset(original_dataset, induction == "Ic" & maintenance == "Mc"))
  n_ItMc <- nrow(subset(original_dataset, induction == "It" & maintenance == "Mc"))
  n_IcMt <- nrow(subset(original_dataset, induction == "Ic" & maintenance == "Mt"))
  n_ItMt <- nrow(subset(original_dataset, induction == "It" & maintenance == "Mt"))
  t_judge <- original_dataset$t_judge[1]

  pseudo_dataset <- generate_scenario(
    sim_num = num,
    t_judge = t_judge,
    censor_rate = censor_rate,
    IcMc_size = n_IcMc,
    ItMc_size = n_ItMc,
    IcMt_size = n_IcMt,
    ItMt_size = n_ItMt,
    theta_Ic = estimators$theta_Ic,
    theta_It = estimators$theta_It,
    lambda_Ic_r = 0.0,
    lambda_Ic_nr = estimators$lambda_Ic_nr,
    lambda_It_r = 0.0,
    lambda_It_nr = estimators$lambda_It_nr,
    lambda_IcMc_r = estimators$lambda_IcMc_r,
    lambda_IcMc_nr = estimators$lambda_IcMc_nr,
    lambda_ItMc_r = estimators$lambda_ItMc_r,
    lambda_ItMc_nr = estimators$lambda_ItMc_nr,
    lambda_IcMt_r = estimators$lambda_IcMt_r,
    lambda_IcMt_nr = estimators$lambda_IcMt_nr,
    lambda_ItMt_r = estimators$lambda_ItMt_r,
    lambda_ItMt_nr = estimators$lambda_ItMt_nr,
    seed = seed
  )
  pseudo_dataset
}

fisher_information_monte_carlo <- function (
  pseudo_datasets,
  estimators
) {
  # second derivative of the lower bound of loglikelihood (matrix)
  qp2_mat <- diag(c(
    estimators$theta_Ic_info,
    estimators$theta_It_info,
    estimators$lambda_Ic_nr_info,
    estimators$lambda_It_nr_info,
    estimators$lambda_IcMc_r_info,
    estimators$lambda_IcMc_nr_info,
    estimators$lambda_ItMc_r_info,
    estimators$lambda_ItMc_nr_info,
    estimators$lambda_IcMt_r_info,
    estimators$lambda_IcMt_nr_info,
    estimators$lambda_ItMt_r_info,
    estimators$lambda_ItMt_nr_info
  ))

  # first derivative of the complete loglikelihood (vector)
  ids <- unique(pseudo_datasets$id)
  grad_func <- function (target_id) {
    dataset <- pseudo_datasets[pseudo_datasets$id == target_id,]
    loglikelihood_grad_vector(dataset, estimators)
  }
  grad_vectors <- lapply(ids, grad_func)

  # approximate variance of the first derivative of observative loglikelihood
  grad_prod <- lapply(grad_vectors, function (l) {l %*% t(l)})
  grad_prod_sum <- Reduce(`+`, grad_prod)

  grad_sum <- Reduce(`+`, grad_vectors)
  grad_sum_prod <- grad_sum %*% t(grad_sum)

  # approximate Fisher's information matrix
  approx_fisher_info <- (
    qp2_mat
    - grad_prod_sum / length(ids)
    + grad_sum_prod / (length(ids)^2)
  )
  approx_fisher_info
}

