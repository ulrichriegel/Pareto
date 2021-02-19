library(Pareto)
context("test functions Pareto")

test_that("Pareto_Layer_Mean", {
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2), 1600)
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 1000), 400)
  expect_equal(Pareto_Layer_Mean(8000, 2000, 2, t = 5000, truncation = 10000), 4666.66666666667)
  expect_equal(Pareto_Layer_Mean(2000, 1000, 0, truncation = 5000, t = 500), 835.16520831955449)
  expect_equal(Pareto_Layer_Mean((1:3) * 8000, (1:3) * 2000, 2), (1:3) * 1600)
  expect_equal(Pareto_Layer_Mean((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 1000), (1:3) * 400)
  expect_equal(Pareto_Layer_Mean((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 1000, truncation = NULL), (1:3) * 400)
  expect_equal(Pareto_Layer_Mean((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 5000, truncation = (1:3) * 10000), (1:3) * 4666.66666666667)
})

test_that("Pareto_Layer_SM", {
  expect_equal(Pareto_Layer_SM(8000, 2000, 2), 6475503.2994728)
  expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 1000), 1618875.8248682)
  expect_equal(Pareto_Layer_SM(8000, 2000, 2, t = 5000, truncation = 10000), 23543145.370663)
  expect_equal(Pareto_Layer_SM(2000, 1000,0, truncation = 5000, t = 1500), 2584318.901925616)
  expect_equal(Pareto_Layer_SM((1:3) * 8000, (1:3) * 2000, 2), (1:3)^2 * 6475503.2994728)
  expect_equal(Pareto_Layer_SM((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 1000), (1:3)^2 * 1618875.8248682)
  expect_equal(Pareto_Layer_SM((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 5000, truncation = (1:3) * 10000), (1:3)^2 * 23543145.370663)
})

test_that("Pareto_Layer_Var", {
  expect_equal(Pareto_Layer_Var(8000, 2000, 2), 3915503.2994728)
  expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
  expect_equal(Pareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)
  expect_equal(Pareto_Layer_Var(2000,1000,0, truncation = 5000, t = 500), 667015.32799764269)
  expect_equal(Pareto_Layer_Var((1:3) * 8000, (1:3) * 2000, 2), (1:3)^2 * 3915503.2994728)
  expect_equal(Pareto_Layer_Var((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 1000), (1:3)^2 * 1458875.8248682)
  expect_equal(Pareto_Layer_Var((1:3) * 8000, (1:3) * 2000, 2, t = (1:3) * 5000, truncation = (1:3) * 10000), (1:3)^2 * 1765367.59288524)

})

test_that("pPareto", {
  expect_equal(pPareto(c(1:3) * 2000, t = (1:3) * 1000, alpha = 2, truncation = NULL), rep(0.75, 3))
  expect_equal(pPareto(c(1:3) * 2000, t = (1:3) * 1000, alpha = 2, truncation = (1:3) * 10000), rep(0.75757575757575757, 3))
})

test_that("dPareto", {
  expect_equal(dPareto(c(1:3) * 2000, t = (1:3) * 1000, alpha = 2, truncation = NULL), c(2.5000000000000001e-04, 1.2500000000000000e-04, 8.3333333333333331e-05))
})

test_that("qPareto", {
  expect_equal(qPareto((1:3) * 0.3, (1:3) * 1000, 1:3), c(1428.5714285714287, 3162.2776601683790, 6463.3040700956490))
  expect_equal(qPareto((1:3) * 0.3, (1:3) * 1000, 1:3, truncation = 10000), c(1369.8630136986301, 3071.4755841697556, 6011.2419963076682))
  expect_equal(qPareto((1:3) * 0.3, (1:3) * 1000, 1:3, truncation = (1:3) * 10000), c(1369.8630136986301, 3138.8241028717225, 6444.0296890428708))
})



################################
# test rPareto with truncation #
################################

test_that("rPareto with truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  truncation <- 20000
  alpha <- 1.5

  expect_equal(Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t, truncation = truncation), 700.14314962210699)

  losses <- rPareto(NumberOfSimulations, t, alpha, truncation = truncation)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t, truncation = truncation), 2)

  expect_equal(ratio, 1)
})

############################################
# test rPareto without truncation #
############################################

test_that("rPareto without truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- 1000
  Cover <- 8000
  AttPoint <- 2000
  alpha <- 1.5

  expect_equal(Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t), 781.75803033941929)

  losses <- rPareto(NumberOfSimulations, t, alpha)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / Pareto_Layer_Mean(Cover, AttPoint,  alpha, t = t), 2)

  expect_equal(ratio, 1)
})


test_that("Pareto_Extrapolation", {
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2), 0.5)
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2, ExpLoss_1 = 1000), 500)
  expect_equal(Pareto_Extrapolation(1000, 1000, 2000, 2000, alpha = 2, truncation = 3000, ExpLoss_1 = 1000), 142.85714285714292)
  expect_equal(Pareto_Extrapolation(Inf, 1000, Inf, 2000, alpha = 0, truncation = 3000, ExpLoss_1 = 100), 20.975411735345432)
  expect_equal(Pareto_Extrapolation(c(1000, 1000, 1000, Inf), 1000, c(2000, 2000, 2000, Inf), 2000, alpha = c(2, 2, 2, 0), truncation = c(Inf, Inf, 3000, 3000), ExpLoss_1 = c(1, 1000, 1000, 100)), c(0.5, 500, 142.85714285714292, 20.975411735345432))
  expect_equal(Pareto_Extrapolation(c(1000, 1000), 1000, 2000, 2000, alpha = 2, truncation = NULL, ExpLoss_1 = c(1, 1000)), c(0.5, 500))
})

test_that("Pareto_Find_Alpha_btw_FQ_Layer", {
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 2000, 2000, 100), 2.9330042139247037)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(2000, 0.25, 1000, 1000, 500), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(1000, 1, 1000, 1000, 500, truncation = 5000), 1.8363401129702193)
  expect_equal(Pareto_Find_Alpha_btw_FQ_Layer(c(2000, 1000), c(0.25, 1), 1000, 1000, 500, truncation = c(Inf, 5000)), c(2, 1.8363401129702193))
})

test_that("Pareto_Find_Alpha_btw_Layers", {
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 100), 1)
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 50), 2)
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 50, truncation = 5000), 1.3871313763147217)
  expect_equal(Pareto_Find_Alpha_btw_Layers(1000, 1000, 100, 2000, 2000, 50, truncation = c(Inf, 5000)), c(2, 1.3871313763147217))
  expect_equal(Pareto_Find_Alpha_btw_Layers(c(1000, 1000), 1000, 100, 2000, 2000, c(100, 50), truncation = NULL), c(1, 2))
})

test_that("Pareto_Find_Alpha_btw_FQs", {
  expect_equal(Pareto_Find_Alpha_btw_FQs(1000, 1, 2000, 0.5), 1)
  expect_equal(Pareto_Find_Alpha_btw_FQs(2000, 0.25, 1000, 1), 2)
  expect_equal(Pareto_Find_Alpha_btw_FQs(2000, 0.25, 1000, 1, truncation = 4000), 1.5849625007211574)
  expect_equal(Pareto_Find_Alpha_btw_FQs(c(1000, 2000, 2000), c(1, 0.25, 0.25), c(2000, 1000, 1000), c(0.5, 1, 1), truncation = c(Inf, Inf, 4000)), c(1, 2, 1.5849625007211574))
})


##################################
# test Pareto_ML_Estimator_alpha #
##################################

test_that("Pareto_ML_Estimator_alpha", {
  losses <- c(1622.49986584698, 1025.1735923535, 1142.67198754259, 1598.2131674777, 1369.79742768744, 1006.5249344124,
              2019.3663238659, 1007.2758879241, 1377.79293040511, 1605.21438984656)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000), 4), 3.4060)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, truncation = 3000), 4), 2.9601)
  losses2 <- c(losses, losses[1:2])
  w <- rep(1, length(losses))
  w[1:2] <- 2
  expect_equal(Pareto_ML_Estimator_Alpha(losses, 1000, weights = w), Pareto_ML_Estimator_Alpha(losses2, 1000))
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, weights = w, truncation = 3000), 5), round(Pareto_ML_Estimator_Alpha(losses2, 1000, truncation = 3000), 5))

  reporting_thresholds <- c(1000, 1000, 1000, 1200, 1200, 1000, 1500, 1000, 1000, 1000)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds), 5), 4.61698)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds, weights = w, truncation = 3000), 5), 4.20722)
  censored <- rep(FALSE, length(losses))
  censored[1:2] <- TRUE
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds, is.censored = censored, weights = w, truncation = 3000), 5), 2.47743)
  expect_equal(round(Pareto_ML_Estimator_Alpha(losses, 1000, reporting_thresholds = reporting_thresholds, is.censored = censored), 5), 3.69358)

  t <- 100
  reporting_thresholds <- rPareto(100, 100, 2)
  alpha <- 2
  losses <- rPareto(100, reporting_thresholds, alpha)
  w <- rep(1, 100)
  w[1:2] <- 10
  losses2 <- c(losses, rep(losses[1:2], 9))
  reporting_thresholds2 <- c(reporting_thresholds, rep(reporting_thresholds[1:2], 9))
  expect_equal(Pareto_ML_Estimator_Alpha(losses, t, reporting_thresholds = reporting_thresholds, weights = w), Pareto_ML_Estimator_Alpha(losses2, t, reporting_thresholds = reporting_thresholds2))
  truncation <- 2 * max(losses)
  expect_equal(Pareto_ML_Estimator_Alpha(losses, t, reporting_thresholds = reporting_thresholds, weights = w, truncation = truncation), Pareto_ML_Estimator_Alpha(losses2, t, reporting_thresholds = reporting_thresholds2, truncation = truncation))
})


###########################
# test Local_Pareto_Alpha #
###########################

test_that("Local_Pareto_Alpha", {
  x <- 1:10 * 1e6
  expect_equal(Local_Pareto_Alpha(x, "norm", mean = 5e6, sd = 2e6), c(0.027623931339495, 0.138789750458851, 0.431399956408768, 1.01832086767407, 1.99471140200716, 3.42323331110419, 5.33797346656343, 7.75470866649017, 10.6794698977028, 14.1137239883195))
  expect_equal(Local_Pareto_Alpha(x, "lnorm", meanlog = 0, sdlog = 4), c(0.92698000118132, 0.967875348611523, 0.991874809629857, 1.00893510366812, 1.02218585767522, 1.03302367422092, 1.04219455399385, 1.05014426920658, 1.0571606056358, 1.0634402085453))
  expect_equal(Local_Pareto_Alpha(x, "gamma", shape = 5, rate = 1e-6), c(0.0153846153846154, 0.19047619047619, 0.618320610687023, 1.24271844660194, 1.99171446781389, 2.81739130434783, 3.69141225565561, 4.59708193041527, 5.52427729441482, 6.46663217796171))
  expect_equal(Local_Pareto_Alpha(x, "weibull", shape = 0.5, scale = 1e6), c(0.5, 0.707106781186547, 0.866025403784439, 1, 1.11803398874989, 1.22474487139159, 1.3228756555323, 1.4142135623731, 1.5, 1.58113883008419))
  expect_equal(Local_Pareto_Alpha(x, "exp", rate = 1e-6), c(1, 2, 3, 3.99999999999999, 4.99999999999999, 5.99999999999996, 7.00000000000033, 7.99999999999989, 9.0000000000022, 9.9999999999987))
  expect_equal(Local_Pareto_Alpha(x, "Pareto", t = 1e6, alpha = 1, truncation = 20e6), c(1.05263157894737, 1.11111111111111, 1.17647058823529, 1.25, 1.33333333333333, 1.42857142857143, 1.53846153846154, 1.66666666666667, 1.81818181818182, 2))
  expect_equal(Local_Pareto_Alpha(x, "GenPareto", t = 1e6, alpha_ini = 1, alpha_tail = 2), c(1, 1.33333333333333, 1.5, 1.6, 1.66666666666667, 1.71428571428571, 1.75, 1.77777777777778, 1.8, 1.81818181818182))
  expect_equal(Local_Pareto_Alpha(x, "PiecewisePareto", t = c(1e6, 3e6, 5e6), alpha = c(1, 2, 3), truncation = 20e6, truncation_type = "wd"), c(1.00187852222918, 1.00376411543287, 2.01131363922062, 2.02020202020202, 3.04761904761905, 3.08324768756423, 3.13438683557529, 3.2051282051282, 3.30078393618485, 3.42857142857144))

})

