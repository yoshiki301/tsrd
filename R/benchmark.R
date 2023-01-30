nonpara_benchmark <- function (
  dataset
) {
  X_Ic <- subset(dataset, induction == "Ic", select = c(time, cens, t_judge, is_induction, maintenance))
  X_It <- subset(dataset, induction == "It", select = c(time, cens, t_judge, is_induction, maintenance))
  t_judge <- dataset$t_judge[1]

  X_Ic_before <- X_Ic
  X_Ic_before$cens <- ifelse(X_Ic$is_induction == 0, 0, X_Ic$cens)
  X_Ic_after <- subset(X_Ic, is_induction == 0)
  X_Ic_after$is_treatment <- ifelse(X_Ic_after$maintenance == "Mt", 1, 0)
  X_It_before <- X_It
  X_It_before$cens <- ifelse(X_It$is_induction == 0, 0, X_It$cens)
  X_It_after <- subset(X_It, is_induction == 0)
  X_It_after$is_treatment <- ifelse(X_It_after$maintenance == "Mt", 1, 0)

  km_fit_Ic <- survival::survfit(survival::Surv(time, cens) ~ 1, data = X_Ic_before)
  km_fit_It <- survival::survfit(survival::Surv(time, cens) ~ 1, data = X_It_before)
  cox_fit_Ic <- survival::coxph(survival::Surv(time - t_judge, cens) ~ is_treatment, data = X_Ic_after)
  cox_fit_It <- survival::coxph(survival::Surv(time - t_judge, cens) ~ is_treatment, data = X_It_after)
  rmst_Ic <- rmst_km(km_fit_Ic, t_judge)
  rmst_It <- rmst_km(km_fit_It, t_judge)
  hr_Ic <- exp(coef(cox_fit_Ic))[[1]]
  hr_It <- exp(coef(cox_fit_It))[[1]]

  list(
    rmst_Ic = rmst_Ic,
    rmst_It = rmst_It,
    hr_Ic = hr_Ic,
    hr_It = hr_It
  )
}

rmst_km <- function (
  surv_fit,
  t_judge
) {
  target_index <- which(surv_fit$time <= t_judge)
  target_time <- c(0, surv_fit$time[target_index], t_judge)
  dt <- diff(target_time)
  target_surv <- c(1.0, surv_fit$surv[target_index])
  rmst_hat <- sum(dt * target_surv)
  rmst_hat
}

para_benchmark <- function (
  dataset
) {
  X_Ic <- subset(dataset, induction == "Ic", select = c(time, cens, t_judge, is_induction, maintenance))
  X_It <- subset(dataset, induction == "It", select = c(time, cens, t_judge, is_induction, maintenance))
  t_judge <- dataset$t_judge[1]

  X_Ic_before <- X_Ic
  X_Ic_before$cens <- ifelse(X_Ic$is_induction == 0, 0, X_Ic$cens)
  X_IcMc_after <- subset(X_Ic, is_induction == 0 & maintenance == "Mc")
  X_IcMt_after <- subset(X_Ic, is_induction == 0 & maintenance == "Mt")
  X_It_before <- X_It
  X_It_before$cens <- ifelse(X_It$is_induction == 0, 0, X_It$cens)
  X_ItMc_after <- subset(X_It, is_induction == 0 & maintenance == "Mc")
  X_ItMt_after <- subset(X_It, is_induction == 0 & maintenance == "Mt")

  reg_fit_Ic <- survival::survreg(survival::Surv(time, cens) ~ 1, data = X_Ic_before, dist = "exponential")
  reg_fit_It <- survival::survreg(survival::Surv(time, cens) ~ 1, data = X_It_before, dist = "exponential")
  lambda_Ic <- (1 / reg_fit_Ic$icoef)[[1]]
  lambda_It <- (1 / reg_fit_It$icoef)[[1]]
  rmst_Ic <- rmst_parametric(lambda_Ic, t_judge)
  rmst_It <- rmst_parametric(lambda_It, t_judge)

  reg_fit_IcMc <- survival::survreg(survival::Surv(time - t_judge, cens) ~ 1, data = X_IcMc_after, dist = "exponential")
  reg_fit_IcMt <- survival::survreg(survival::Surv(time - t_judge, cens) ~ 1, data = X_IcMt_after, dist = "exponential")
  reg_fit_ItMc <- survival::survreg(survival::Surv(time - t_judge, cens) ~ 1, data = X_ItMc_after, dist = "exponential")
  reg_fit_ItMt <- survival::survreg(survival::Surv(time - t_judge, cens) ~ 1, data = X_ItMt_after, dist = "exponential")
  lambda_IcMc <- (1 / reg_fit_IcMc$icoef)[[1]]
  lambda_IcMt <- (1 / reg_fit_IcMt$icoef)[[1]]
  lambda_ItMc <- (1 / reg_fit_ItMc$icoef)[[1]]
  lambda_ItMt <- (1 / reg_fit_ItMt$icoef)[[1]]
  hr_Ic <- lambda_IcMt / lambda_IcMc
  hr_It <- lambda_ItMt / lambda_ItMc

  list(
    rmst_Ic = rmst_Ic,
    rmst_It = rmst_It,
    hr_Ic = hr_Ic,
    hr_It = hr_It
  )
}

rmst_parametric <- function (
  lambda,
  t_judge
) {
  stats::pexp(t_judge, lambda) / lambda
}

benchmark <- function (
  dataset
) {
  nonpara <- nonpara_benchmark(dataset)
  para <- para_benchmark(dataset)
  list(
    rmst_Ic_km = nonpara$rmst_Ic,
    rmst_It_km = nonpara$rmst_It,
    hr_Ic_cox = nonpara$hr_Ic,
    hr_It_cox = nonpara$hr_It,
    rmst_Ic_exp = para$rmst_Ic,
    rmst_It_exp = para$rmst_It,
    hr_Ic_exp = para$hr_Ic,
    hr_It_exp = para$hr_It
  )
}

add_benchmark_column <- function (
  datasets,
  treatment_effects
) {
  ids <- unique(datasets$id)
  add_func <- function (target_id) {
    dataset <- datasets[datasets$id == target_id,]
    benchmark_col <- benchmark(dataset)
    treatment_effect <- treatment_effects[treatment_effects$id == target_id,]
    cbind(treatment_effect, benchmark_col)
  }

  add_list <- lapply(ids, add_func)
  Reduce(rbind, add_list)
}
