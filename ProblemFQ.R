load("BCUNL.dat")
load("AP.dat")

Result <- PiecewisePareto_Match_Layer_Losses(attachment_points, BC_unl, Unlimited_Layers = T)
Result
PiecewisePareto_CDF(1e9, Result$t, Result$alpha)
min(Result$alpha)
PiecewisePareto_Layer_Mean(Inf, 0, Result$t, Result$alpha)
