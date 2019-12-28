library(Pareto)
context("test functions")

test_that("PiecewisePareto_Layer_Mean", {
  expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, 2, t = 1000), 400)
  expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, 2, t = 5000, truncation = 10000), 4666.66666666667)
})

test_that("PiecewisePareto_Layer_SM", {
  expect_equal(PiecewisePareto_Layer_SM(8000, 2000, 2, t = 1000), 1618875.8248682)
  expect_equal(PiecewisePareto_Layer_SM(8000, 2000, 2, t = 5000, truncation = 10000), 23543145.370663)
})

test_that("PiecewisePareto_Layer_SM", {
  expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
  expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)
})


#####################################################
# test rPiecewisePareto with truncation_type = "wd" #
#####################################################

test_that("rPiecewisePareto with truncation wd", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- c(1000,3000,5000)
  Cover <- 8000
  AttPoint <- 2000
  truncation <- 6000
  alpha <- c(0.7,1.5,2)
  truncation_type ="wd"

  expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 868.729076663466)

  losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 2)

  expect_equal(ratio, 1)
})

#####################################################
# test rPiecewisePareto with truncation_type = "lp" #
#####################################################

test_that("rPiecewisePareto with truncation lp", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- c(1000,3000,5000)
  Cover <- 8000
  AttPoint <- 2000
  truncation <- 6000
  alpha <- c(0.7,1.5,2)
  truncation_type ="lp"

  expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 1255.52081303827)

  losses <- rPiecewisePareto(NumberOfSimulations, t, alpha, truncation = truncation, truncation_type = truncation_type)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha, truncation = truncation, truncation_type = truncation_type), 2)

  expect_equal(ratio, 1)
})

############################################
# test rPiecewisePareto without truncation #
############################################

test_that("rPiecewisePareto without truncation", {
  set.seed(1972)
  NumberOfSimulations <- 1000000

  t <- c(1000,3000,5000)
  Cover <- 8000
  AttPoint <- 2000
  alpha <- c(0.7,1.5,2)

  expect_equal(PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha), 1696.1079667876)

  losses <- rPiecewisePareto(NumberOfSimulations, t, alpha)
  xs_losses <- pmin(Cover, pmax(0, losses - AttPoint))
  mean(xs_losses)

  ratio <- round(mean(xs_losses) / PiecewisePareto_Layer_Mean(Cover, AttPoint,  t, alpha), 2)

  expect_equal(ratio, 1)
})


###########################################
# test PiecewisePareto_ML_Estimator_alpha #
###########################################

test_that("Pareto_ML_Estimator_alpha", {
  losses <- c(2021.37427270114, 1144.4154953111, 1263.31435215882, 2159.17623012884, 1597.37248664761, 1272.28453257,
              2914.10585100956, 1304.18266162408, 1120.50113120055, 2127.70698437685, 2060.81683374026, 2039.7843971136,
              1888.54098901941, 2124.60113577822, 1861.90023645773, 1061.42587835879, 1096.87195504782, 2741.00486150546,
              3599.51085744179, 1944.00791988792, 4801.21672397748, 2011.17712322692, 1614.38016727125, 1366.00934974079,
              2740.17455761229, 3696.80509791037, 1519.32575139179, 2673.72771069715, 1377.97530890359, 2747.57267315656,
              2962.51959921929, 2970.06262726659, 3849.91766409146, 1011.07955688529, 2142.35280378138, 2363.88639730512,
              1585.05201649176, 1519.87964078085, 4152.85247579578, 1145.85174532018, 1045.70673772088, 2035.21218839728,
              2265.69458584306, 1259.17343777875, 1276.43540240796, 1156.59174646013, 2213.19952527109, 2090.8925597549,
              1599.32662288134, 1748.10607173994, 1110.83766088527, 4216.9888129572, 1987.10180900716, 2848.75495497778,
              2459.33207487438, 1057.30210620244, 2613.23427141216, 4267.61547888255, 2110.3988745262, 1608.79852597768,
              3299.23731214781, 1598.18841764305, 1867.71915473237, 1025.8298666397, 1714.23234124254, 2344.71364141762,
              1415.68319510445, 1529.97435270215, 1236.57055448752, 2466.0269323041, 1091.02769697244, 2408.02293167307,
              3106.53753930817, 2461.21276690245, 3262.03972168424, 4127.01858195009, 1181.45647693256, 2915.05440240728,
              2819.40842718265, 1079.96118940952, 3652.61964913811, 3731.17802733447, 4689.68469988139, 3430.45722285605,
              1700.45433736218, 2269.02625009652, 3545.91823730921, 3112.72960104325, 2757.85209289505, 1073.92628882837,
              3697.84111140716, 2567.88718258095, 2624.17812649749, 2430.48340413282, 2399.39130252127, 1194.87342560877,
              1733.25071770226, 2177.58393770677, 2313.84456911248, 2561.30080556934)
  expect_equal(round(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000)), 4), c(0.8126, 2.6572, 4.4220))
  expect_equal(round(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000), truncation = 5000, truncation_type = "lp"), 4), c(0.8126, 2.6572, 1.3569))
  expect_equal(round(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000), truncation = 5000, truncation_type = "wd"), 4), c(0.5965, 1.5063, 0.9891))
})


