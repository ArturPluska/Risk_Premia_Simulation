# Transfering from Matrix to List

MatrixOrFrameToList <- function(dataSeries) {
  if(is.matrix(dataSeries) || is.data.frame(dataSeries)) {
    matrixToList <- list()
    for(i in 1:ncol(dataSeries)) {
      matrixToList <- c(matrixToList, list(dataSeries[, i, drop = T]))
    }
    names(matrixToList) <- colnames(dataSeries)
  } else {
    stop('Wrong "dataSeries" input! It should be a matrix.')
  }
  return(matrixToList)
}

# Rate of Return

RateReturn <- function(dataSeries, time = 'continuous') {
  if((is.vector(dataSeries) || is.matrix(dataSeries) || is.data.frame(dataSeries))) {
    if((time == 'continuous' || time == 'digital')) {
      if(time == 'continuous') {
        if(is.vector(dataSeries)) {
          changeRate <- log(dataSeries[-1] / dataSeries[-length(dataSeries)])
        } else if(is.matrix(dataSeries)) {
          changeRate <- matrix(NA, nrow = dim(dataSeries)[1] - 1, ncol = dim(dataSeries)[2])
          colnames(changeRate) <- colnames(dataSeries)
          for(i in 1:dim(dataSeries)[2]) {
            for(j in 1:dim(dataSeries)[1] - 1) {
              changeRate[j, i] <- log(dataSeries[j + 1, i]/ dataSeries[j , i])
            }
          }
        } else if(is.data.frame(dataSeries)) {
          changeRate <- matrix(NA, nrow = dim(dataSeries)[1] - 1, ncol = dim(dataSeries)[2])
          colnames(changeRate) <- colnames(dataSeries)
          changeRate <- as.data.frame(changeRate)
          for(i in 1:dim(dataSeries)[2]) {
            for(j in 1:dim(dataSeries)[1] - 1) {
              changeRate[j, i] <- log(dataSeries[j + 1, i]/ dataSeries[j , i])
            }
          }
        }
      } else if(time == 'digital') {
        if(is.vector(dataSeries)) {
          changeRate <- dataSeries[-1] / dataSeries[-length(dataSeries)]
        } else if(is.matrix(dataSeries)) {
          changeRate <- matrix(NA, nrow = dim(dataSeries)[1] - 1, ncol = dim(dataSeries)[2])
          colnames(changeRate) <- colnames(dataSeries)
          for(i in 1:dim(dataSeries)[2]) {
            for(j in 1:dim(dataSeries)[1] - 1) {
              changeRate[j, i] <- dataSeries[j + 1, i]/ dataSeries[j , i]
            }
          }
        } else if(is.data.frame()) {
          changeRate <- matrix(NA, nrow = dim(dataSeries)[1] - 1, ncol = dim(dataSeries)[2])
          colnames(changeRate) <- colnames(dataSeries)
          changeRate <- as.data.frame(changeRate)
          for(i in 1:dim(dataSeries)[2]) {
            for(j in 1:dim(dataSeries)[1] - 1) {
              changeRate[j, i] <- dataSeries[j + 1, i]/ dataSeries[j , i]
            }
          }
        }
      }  
    } else {
      stop('Wrong "time" input! Choose between "continuous" and "digital".')
    }
  } else {
    stop('Wrong "dataSeries" input! It should be as vector or matrix.')
  }
  return(changeRate)
}

# Geometric Brownian Motion (GBM) with Gauss

GeometricBrownianMotionGS <- function(startPrice, mu, sigma, steps, horizon) {
  dt <- horizon / steps
  gbmgs <- startPrice * exp((mu - sigma^2 / 2) * cumsum(c(0, rep(dt, steps))) + c(0, sigma * cumsum(sqrt(dt) * rnorm(steps, 0, 1))))
  return(gbmgs)
}

# Vasicek with Gauss

VasicekGS <- function(startRate, theta, mu, sigma, steps, horizon) {
  dt <- horizon / steps
  rate <- startRate
  vgs <- c(startRate)
  for(i in 1:steps) {
    dr <- theta * (mu - rate) * dt + sigma * rnorm(1, 0, 1)
    rate <- rate + dr
    vgs <- c(vgs, rate)
  }
  return(vgs)
}

# Cox Ingersoll Ross Model (CIR) with Gauss

CoxIngersollRossGS <- function(startRate, theta, mu, sigma, steps, horizon) {
  dt <- horizon / steps
  rate <- startRate
  cirgs <- c(startRate)
  for(i in 1:steps) {
    dr <- theta * (mu - rate) * dt + sigma * sqrt(rate) * rnorm(1, 0, 1)
    rate <- ifelse(rate + dr < 0.001, 0.001, rate + dr)
    cirgs <- c(cirgs, rate)
  }
  return(cirgs)
}

# Correlated GBM and Vasicek with Gauss

CorGBMwVAS <- function(startStock, stockMu, stockSigma, startRate, rateTheta, rateMu, rateSigma, corrCoef, steps, horizon) {
  dt <- horizon / steps
  stoch <- matrix(rnorm(2 * steps, 0, 1), nrow = steps, ncol = 2) %*% chol(matrix(c(1, corrCoef, corrCoef, 1), 2, 2))
  stock <- startStock
  gbm <- c(startStock)
  for(i in 1:steps) {
    ret <- exp((stockMu - stockSigma^2 / 2) * dt + stockSigma * sqrt(dt) * stoch[i, 1])
    stock <- stock * ret
    gbm <- c(gbm, stock)
  }
  rate <- startRate
  vas <- c(startRate)
  for(i in 1:steps) {
    drate <- rateTheta * (rateMu - rate) * dt + rateSigma * stoch[i, 2]
    rate <- rate + drate
    vas <- c(vas, rate)
  }
  cgbmv <- cbind(gbm, vas)
  colnames(cgbmv) <- c('Equity', 'FI')
  rownames(cgbmv) <- c(1:(steps+1))
  return(as.data.frame(cgbmv, row.names = NULL))
}

# Correlated GBM and CIR with Gauss

CorGBMwCIR <- function(startStock, stockMu, stockSigma, startRate, rateTheta, rateMu, rateSigma, corrCoef, steps, horizon) {
  dt <- horizon / steps
  stoch <- matrix(rnorm(2 * steps, 0, 1), nrow = steps, ncol = 2) %*% chol(matrix(c(1, corrCoef, corrCoef, 1), 2, 2))
  stock <- startStock
  gbm <- c(startStock)
  for(i in 1:steps) {
    ret <- exp((stockMu - stockSigma^2 / 2) * dt + stockSigma * sqrt(dt) * stoch[i, 1])
    stock <- stock * ret
    gbm <- c(gbm, stock)
  }
  rate <- startRate
  cir <- c(startRate)
  for(i in 1:steps) {
    drate <- rateTheta * (rateMu - rate) * dt + rateSigma * sqrt(dt) * stoch[i, 2]
    rate <- ifelse(rate + drate < 0.001, 0.001, rate + drate)
    cir <- c(cir, rate)
  }
  cgbmcir <- cbind(gbm, cir)
  colnames(cgbmcir) <- c('Equity', 'FI')
  rownames(cgbmcir) <- c(1:(steps+1))
  return(as.data.frame(cgbmcir, row.names = NULL))
}

# Monte Carlo Simulation GBM

MonteCarloSimulationEquity <- function(simulationModel, parametersVector, simulationNumber) {
  monteCarloScenarios <- sapply(X = rep(parametersVector$startPrice, simulationNumber),
                                FUN = simulationModel,
                                mu = parametersVector$mu,
                                sigma = parametersVector$sigma,
                                steps = parametersVector$steps, 
                                horizon = parametersVector$horizon)
  return(as.data.frame(monteCarloScenarios))
}

# Monte Carlo Simulation Vasicek or CIR

MonteCarloSimulationFixedIncome <- function(simulationModel, parametersVector, simulationNumber) {
  monteCarloScenarios <- sapply(X = rep(parametersVector$startRate, simulationNumber),
                                FUN = simulationModel,
                                theta = parametersVector$theta,
                                mu = parametersVector$mu,
                                sigma = parametersVector$sigma,
                                steps = parametersVector$steps, 
                                horizon = parametersVector$horizon)
  return(as.data.frame(monteCarloScenarios))
}

# Monte Carlo Simulation BM with Vasicek or CIR

MonteCarloSimulationEquityFixedIncome <- function(simulationModel, parametersVector, simulationNumber) {
  tmp <- sapply(X = rep(parametersVector$steps, simulationNumber),
                                FUN = simulationModel,
                                startStock = parametersVector$startStock,
                                stockMu = parametersVector$stockMu,
                                stockSigma = parametersVector$stockSigma,
                                rateTheta = parametersVector$rateTheta,
                                startRate = parametersVector$startRate,
                                rateMu = parametersVector$rateMu,
                                rateSigma = parametersVector$rateSigma,
                                corrCoef = parametersVector$corrCoef,
                                horizon = parametersVector$horizon)
  monteCarloScenarios <- list()
  nameSeries <- c('Equity', 'FI')
  for(nm in nameSeries) {
    monteCarloScenarios[[nm]] <- as.data.frame(sapply(X = c(1:simulationNumber), function(x) {as.data.frame(tmp)[[x]][nm]}), col.names = paste('S', c(1:simulationNumber), sep  = '')) 
  }
  return(monteCarloScenarios)
}

# Zero Coupon Bond Price under Vasicek model

ZeroCouponBondPriceVasicek <- function(rate, theta, mu, sigma, tenor) {
  a_par <- (1 / theta) * (1 - exp(- tenor * theta))
  b_par <- (mu - (sigma^2 / (2 * theta ^ 2))) * (a_par - tenor) - ((sigma^2 / (4 * theta)) * a_par^2)
  return(exp(b_par - a_par * rate))
}

# Zero Coupon Bond (ZCB) Price under CIR model

#tbc

# Yield To Maturity ZCB
YieldToMaturityZCB <- function(price, tenor) {
  return(- log(price) / tenor)
}

# Vasicek Calibration by OLS

VasicekCalibrationOLS <- function(alfaCons, betaTrend, residuaOLS, steps, horizon) {
  dt <- horizon / steps
  theta <- ((1 - betaTrend) / dt)
  mu <- (betaTrend / (theta * dt))
  sigma <- sqrt(var(residuaOLS) / dt)
  modelParameters <- c(theta, mu, sigma)
  names(modelParameters) <- c('theta', 'mu', 'sigma')
  return(modelParameters)
}

# CIR Calibration by OLS

#tbc