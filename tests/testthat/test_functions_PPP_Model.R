library(Pareto)
context("test functions Match_Layer_Losses and PPP_Model")

test_that("PiecewisePareto_Match_Layer_Losses", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
})


test_that("PiecewisePareto_Match_Layer_Losses truncated lp", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
})

test_that("PiecewisePareto_Match_Layer_Losses truncated wd", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
})


test_that("PiecewisePareto_Match_Layer_Losses truncated lp & frequencies", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
  expect_equal(PPP_Model_Excess_Frequency(c(1000, 2000, 5000), Fit), c(1.1, 0.95, 0.5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, FQ_at_lowest_AttPt = 1.2, FQ_at_highest_AttPt = 0.4, truncation = 10000)
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
  expect_equal(PPP_Model_Excess_Frequency(c(1000, 2000, 5000), Fit), c(1.2, 0.95, 0.4))

})

test_that("PiecewisePareto_Match_Layer_Losses truncated wd & frequencies", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
  expect_equal(PPP_Model_Excess_Frequency(c(1000, 2000, 5000), Fit), c(1.1, 0.95, 0.5))

  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, FQ_at_lowest_AttPt = 1.2, FQ_at_highest_AttPt = 0.4, truncation = 10000, truncation_type = 'wd')
  expect_equal(is.valid.PPP_Model(Fit), TRUE)
  expect_equal(PPP_Model_Exp_Layer_Loss(Cover, AP, Fit), EL)
  expect_equal(PPP_Model_Excess_Frequency(c(1000, 2000, 5000), Fit), c(1.2, 0.95, 0.4))

})

test_that("PiecewisePareto_Match_Layer_Losses only two layers", {
  Fit <- PiecewisePareto_Match_Layer_Losses(c(1000,2000), c(100,100))
  expect_equal(Fit$alpha, 2)
  Fit <- PiecewisePareto_Match_Layer_Losses(c(1000,2000), c(100,100), truncation = 10000, truncation_type = "wd")
  expect_equal(Fit$alpha, 1.4364891695020006)

})

test_that("Panjer distribution", {
  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 1)
  expect_equal(PPP_Model_Exp_Layer_Loss(1, 0, Fit), 10)
  expect_equal(PPP_Model_Layer_Var(1, 0, Fit), 10)

  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 0.5)
  expect_equal(PPP_Model_Exp_Layer_Loss(1, 0, Fit), 10)
  expect_equal(PPP_Model_Layer_Var(1, 0, Fit), 5)

  Fit <- PPP_Model(FQ = 10, t = 1000, alpha = 2, dispersion = 2)
  expect_equal(PPP_Model_Exp_Layer_Loss(1, 0, Fit), 10)
  expect_equal(PPP_Model_Layer_Var(1, 0, Fit), 20)

})


test_that("Simulation, Sd and Var, dispersion = 1", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd')

  set.seed(1972)

  losses <- PPP_Model_Simulate(100000, Fit)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses), -1), round(PPP_Model_Layer_Sd(1000, 2000, Fit), -1))

  expect_equal(round(PPP_Model_Layer_Sd(1000, 2000, Fit), 3), 939.264)
  expect_equal(round(PPP_Model_Layer_Var(1000, 2000, Fit), 3), 882217.474)
})

test_that("Simulation, Sd and Var, dispersion = 0.63", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd', dispersion = 0.63)

  set.seed(1972)

  losses <- PPP_Model_Simulate(100000, Fit)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses), -1), round(PPP_Model_Layer_Sd(1000, 2000, Fit), -1))

  expect_equal(round(PPP_Model_Layer_Sd(1000, 2000, Fit), 3), 780.873)
  expect_equal(round(PPP_Model_Layer_Var(1000, 2000, Fit), 3), 609762.928)
})


test_that("Simulation, Sd and Var, dispersion = 2", {
  AP <- c(1000, 2000, 3000, 4000, 5000)
  EL <- c(1000, 900, 800, 600, 500)
  Cover <- c(diff(AP), Inf)
  FQs <- c(1.1, 0.95, NA, NA, 0.5)
  Fit <- PiecewisePareto_Match_Layer_Losses(AP, EL, Frequencies = FQs, truncation = 10000, truncation_type = 'wd', dispersion = 2)

  set.seed(1972)

  losses <- PPP_Model_Simulate(100000, Fit)
  layer_losses <- matrix(pmin(1000, pmax(losses - 2000, 0)), nrow = 100000)
  annual_losses <- apply(layer_losses, 1, function(x) sum(x, na.rm = T))

  expect_equal(round(sd(annual_losses) - PPP_Model_Layer_Sd(1000, 2000, Fit), -1), 0)

  expect_equal(round(PPP_Model_Layer_Sd(1000, 2000, Fit), 3), 1272.235)
  expect_equal(round(PPP_Model_Layer_Var(1000, 2000, Fit), 3), 1618581.110)
})
