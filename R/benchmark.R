coxph_induction_rmst <- function (
  dataset
) {
  induction_data <- dataset
  induction_data$cens <- ifelse(induction_data$is_induction == 1, induction_data$cens, 0)
  Ic_data <- subset(induction_data, induction == "Ic")
  It_data <- subset(induction_data, induction == "It")
  Ic_fit <- survival::survreg(survival::Surv(time, cens) ~ 1, data = Ic_data, dist = "exponential")
  It_fit <- survival::survreg(survival::Surv(time, cens) ~ 1, data = It_data, dist = "exponential")
  list(Ic = Ic_fit, It = It_fit)
}

calc_parametric_rmst <- function (
  lambda,
  t_judge
) {
  stats::pexp(t_judge, lambda) / lambda
}

coxph_maintenance_hr <- function (
  dataset
) {
  maintenance_data <- subset(dataset, is_induction == 0)
  maintenance_data$is_treatment <- ifelse(maintenance_data$maintenance == "Mt", 1, 0)
  fit <- survival::coxph(survival::Surv(time - t_judge, cens) ~ is_treatment, data = maintenance_data)
  fit
}
