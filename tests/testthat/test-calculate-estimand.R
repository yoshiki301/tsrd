test_that("rmst works", {
  expect_equal( # scalar
    rmst(0.5, 1.0, 6),
    3.498761,
    tolerance = 1e-5
  )
  expect_equal( # vector
    rmst(rep(0.5, 5L), rep(1.0, 5L), rep(6, 5L)),
    rep(3.498761, 5L),
    tolerance = 1e-5
  )
})

test_that("hr works", {
  expect_equal( # scalar
    hr(0.6, 0.8),
    0.75,
    tolerance = 1e-5
  )
  expect_equal( # vector
    hr(rep(0.6, 5L), rep(0.8, 5L)),
    rep(0.75, 5L),
    tolerance = 1e-5
  )
})
