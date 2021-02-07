library(Pareto)
context("test functions PiecewisePareto")

test_that("PiecewisePareto_Layer_Mean", {
  expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, t = 1000, alpha = 2), 400)
  expect_equal(PiecewisePareto_Layer_Mean(8000, 2000, t = 5000, alpha = 2, truncation = 10000), 4666.66666666667)
  expect_equal(PiecewisePareto_Layer_Mean(c(8000, 2000), c(2000, 1000), t = 5000, alpha = 2, truncation = 10000), c(4666.66666666667, 2000))
  expect_equal(PiecewisePareto_Layer_Mean(c(8000, 2000), c(2000, 1000), t = c(1000, 3000, 5000), alpha = c(1, 2, 3), truncation = 10000), c(976.89367953673514, 1098.61228866810916))
})

test_that("PiecewisePareto_Layer_SM", {
  expect_equal(PiecewisePareto_Layer_SM(8000, 2000, t = 1000, alpha = 2), 1618875.8248682)
  expect_equal(PiecewisePareto_Layer_SM(8000, 2000, t = 5000, alpha = 2, truncation = 10000), 23543145.370663)
  expect_equal(PiecewisePareto_Layer_SM(c(1000, 8000), c(1000, 2000), t = c(1000, 3000, 5000), alpha = c(1, 2, 3), truncation = 10000), c(613705.63888010941, 3300236.16730614332))
})

test_that("PiecewisePareto_Layer_Var", {
  expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 1000), 1458875.8248682)
  expect_equal(PiecewisePareto_Layer_Var(8000, 2000, 2, t = 5000, truncation = 10000), 1765367.59288524)
  expect_equal(PiecewisePareto_Layer_Var(c(1000, 8000), c(1000, 2000), t = c(1000, 3000, 5000), alpha = c(1, 2, 3), truncation = 10000), c(133252.62496190792, 2345914.90618732199))
})

test_that("qPiecewisePareto", {
  expect_equal(qPiecewisePareto((1:9) * 0.1, c(1000,2000), c(1,2)), c(1111.1111111111111, 1250.0000000000000,
            1428.5714285714287, 1666.6666666666667, 2000.0000000000000, 2236.0679774997902,
            2581.9888974716114, 3162.2776601683800, 4472.1359549995805))
})

test_that("pPiecewisePareto", {
  expect_equal(pPiecewisePareto(c(1:4) * 1000, t = (1:3) * 1000, alpha = 1:3, truncation = 4000), c(0, 0.5, 0.77777777777777779, 1))
})

test_that("dPiecewisePareto", {
  expect_equal(dPiecewisePareto(c(1:4) * 1000, t = (1:3) * 1000, alpha = 1:3, truncation = 4000), c(0, 0.0005, 0.00022222222222222221, 0))
})

test_that("dPiecewisePareto", {
  expect_equal(qPiecewisePareto(c(1:4) * 0.2, t = (1:3) * 1000, alpha = 1:3, truncation = 4000), c(1250.0000000000000, 1666.6666666666667, 2236.0679774997902, 3060.1459631738917))
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

  w <- rep(1, length(losses))
  w[1:2] <- 2
  w[3] <- 5
  losses2 <- c(losses, losses[1:2], rep(losses[3], 4))

  expect_equal(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000), weights = w), PiecewisePareto_ML_Estimator_Alpha(losses2, c(1000, 2000, 3000)))
  expect_equal(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000), weights = w, truncation = 5000, truncation_type = "lp"), PiecewisePareto_ML_Estimator_Alpha(losses2, c(1000, 2000, 3000), truncation = 5000, truncation_type = "lp"))
  expect_equal(PiecewisePareto_ML_Estimator_Alpha(losses, c(1000, 2000, 3000), weights = w, truncation = 5000, truncation_type = "wd"), PiecewisePareto_ML_Estimator_Alpha(losses2, c(1000, 2000, 3000), truncation = 5000, truncation_type = "wd"))

})


