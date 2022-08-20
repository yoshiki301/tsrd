generate_scenario <- function (
  sim_num = 1000L,
  t_judge = 6,
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
  seed = 42
) {
  set.seed(seed)

  generate_func <- function (sim_id) {
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

    data_IcMc <- data.frame(
      id = sim_id,
      time = t_IcMc,
      cens = rbinom(IcMc_size, 1L, 1.0 - censor_rate),
      t_judge = t_judge,
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
      time = t_ItMc,
      cens = rbinom(ItMc_size, 1L, 1.0 - censor_rate),
      t_judge = t_judge,
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
      time = t_IcMc,
      cens = rbinom(IcMt_size, 1L, 1.0 - censor_rate),
      t_judge = t_judge,
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
      time = t_ItMt,
      cens = rbinom(ItMt_size, 1L, 1.0 - censor_rate),
      t_judge = t_judge,
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
