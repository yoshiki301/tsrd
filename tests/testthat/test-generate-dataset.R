test_that("generate_scenario works", {
  scenario <- generate_scenario(
    sim_num = 10L,
    IcMc_size = 20L,
    IcMt_size = 30L,
    ItMc_size = 40L,
    ItMt_size = 50L
  )

  time <- scenario$time
  expect_equal(sum(is.na(time)) == 0, T)

  expected_count <- data.frame(
    induction = factor(c("Ic", "It", "Ic", "It")),
    maintenance = factor(c("Mc", "Mc", "Mt", "Mt")),
    Freq = c(200L, 400L, 300L, 500L)
  )
  actual_count <- as.data.frame(
    table(scenario[,c("induction", "maintenance")])
  )
  expect_equal(actual_count, expected_count)
})

test_that("adjust_censored_time works", {
  time <- runif(100L, min = 0.0, max = 1.0)
  cens <- c(rep(1, 90L), rep(0, 10L))

  adjusted_time <- adjust_censored_time(time, cens)
  expect_equal(all(adjusted_time[1:90] == time[1:90]), T)
  expect_equal(all(adjusted_time[91:100] <= time[91:100]), T)
})
