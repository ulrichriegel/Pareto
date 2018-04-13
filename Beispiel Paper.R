library(Pareto)

Thresholds <- c(1000,2000,3000,5000,10000)
ExpLosses <- c(200,100,100,50,50)
LowestFQ <- 0.25

Result <- PiecewisePareto_Match_Layer_Losses(Attachment_Points = Thresholds, Expected_Layer_Losses = ExpLosses, FQ_at_lowest_AttPt = LowestFQ, minimize_ratios = T, alpha_max = 100, Use_unlimited_Layer_for_FQ = TRUE)
Result
g <- numeric()
for (i in 1:9) {
  g <- c(g, (1-PiecewisePareto_CDF(Result$t[i], Result$t, Result$alpha)) * Result$FQ)
}



View(matrix(c(Result$t, Result$alpha, g), nrow=3, byrow = T))

AnzahlStuetzstellen <- 1101
UpperLimit <- 12000
Step <- (UpperLimit - Thresholds[1]) / (AnzahlStuetzstellen - 1)

t <- Thresholds[1] + (1:AnzahlStuetzstellen) * Step - Step

FQ <- numeric(AnzahlStuetzstellen)
EL <- numeric(AnzahlStuetzstellen)
for (i in 1:AnzahlStuetzstellen) {
  FQ[i] <- (1 - PiecewisePareto_CDF(t[i],Result$t, Result$alpha)) * Result$FQ
  EL[i] <- PiecewisePareto_Layer_Mean("unl", t[i], Result$t, Result$alpha) * Result$FQ
}

df <- as.data.frame(matrix(c(t, EL, FQ), ncol=3))
View(df)

write.table(df, "Ergebnis.csv", sep = ";")
#FQ[1000] <- PiecewisePareto_CDF(t[1000],Result$t, Result$alpha) * Result$FQ
