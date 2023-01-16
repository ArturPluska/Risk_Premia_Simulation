# Loading data

wig <- read.csv('wig.csv', sep = ',', dec = '.', header = FALSE)
wigReturn <- RateReturn(wig)

log(wig$V1[181]/wig$V1[1]) # stopa zwrotu od poczatku 85%

plot(wig$V1, type = 'l')

wibor <- read.csv('wibor.csv', sep = ',', dec = '.', header = FALSE) / 100
wiborX <- wibor$V1[1:(dim(wibor)[1]-1)]
wiborY <- wibor$V1[2:dim(wibor)[1]]

# Calibration

equity_Mu <- mean(wigReturn$V1) # 0.005
equity_Sigma <- sd(wigReturn$V1) # 0.06

wiborOLS <- lm(wiborY ~ wiborX)
summary(wiborOLS)
vasicekParameters <- VasicekCalibrationOLS(alfaCons = wiborOLS$coefficients[[1]],
                                          betaTrend = wiborOLS$coefficients[[2]], 
                                          residuaOLS = wiborOLS$residuals,
                                          steps = 180,
                                          horizon = 15) # theta = 0.15, mu = 0.035, sigma = 0.006
# Simulations

sim_GBM_VAS <- MonteCarloSimulationEquityFixedIncome(simulationModel = CorGBMwVAS,
                                              parametersVector = list(startStock = 100,
                                                                      stockMu = 0.005,
                                                                      stockSigma = 0.06,
                                                                      startRate = 0.015,
                                                                      rateTheta = 0.15,
                                                                      rateMu = 0.035,
                                                                      rateSigma = 0.006,
                                                                      corrCoef = -0.2,
                                                                      steps = 120,
                                                                      horizon = 10),
                                              simulationNumber = 1000)

# Returns and Yields

rate_mean_MC <- c()
for(i in 2:121) {
  rate_mean_MC <- c(rate_mean_MC, mean(as.numeric(sim_GBM_VAS$FI[i, ])))
}

# Visualisation
p1 <- ggplot(data = data.frame(Time = c(1:length(rate_mean_MC)), MeanRate = rate_mean_MC), aes(x = Time, y = MeanRate)) +
        geom_line(colour = 'red', linewidth = 1, linetype = 'solid') + 
        ggtitle("Evolution of Mean Rate") + 
        xlab("Time") + 
        ylab("Mean Rate")
p1

equity_mean_MC <- c()
for(i in 1:109) {
  equity_mean_MC <- c(equity_mean_MC, mean(as.numeric(log(sim_GBM_VAS$Equity[i + 12, ] / sim_GBM_VAS$Equity[i, ]))))
}

p2 <- ggplot(data = data.frame(Time = c(1:length(equity_mean_MC)), MeanEquity = equity_mean_MC), aes(x = Time, y = MeanEquity)) +
        geom_line(colour = 'blue', linewidth = 1, linetype = 'solid') + ggtitle("Evolution of Equity Return") + 
        xlab("Time") + 
        ylab("Equity Return")
p2

p3 <- ggplot(data = data.frame(Time = c(1:length(equity_mean_MC)), Premia = (equity_mean_MC - rate_mean_MC[12:120])), aes(x = Time, y = Premia)) +
  geom_line(colour = 'violet', linewidth = 1, linetype = 'solid') + ggtitle("Premia") + 
  xlab("Time") + 
  ylab("Premia")
p3