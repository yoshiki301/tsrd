test_that("generate_scenario works", {
  scenario <- generate_scenario(sim_num = 10L)

  time <- scenario$time
  expect_equal(sum(is.na(time)) == 0, T)

  expected_count <- data.frame(
    induction = c("Ic", "It", "Ic", "It"),
    maintenance = c("Mc", "Mc", "Mt", "Mt"),
    Freq = c(10000L, 10000L, 10000L, 10000L)
  )
  actual_count <- as.data.frame(
    table(scenario[,c("induction", "maintenance")])
  )
  expect_equal(actual_count, expected_count)
})
