# simulate WHP with inhomogenous baseline rate
# Baseline rate is mu + B * s(t) - C * sin(t/P)
# where B, C are non-negative and mu > C
# and s(t) is 1 when pi<t<2*pi, 3*pi<t<4*pi, 5*pi<t<6*pi etc

# following block allows the script to be run in parallel on NeSI
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
set.seed(4321 + 10*i)

# Integral of baseline

Mu <- function(t, pars){
  B <- pars[["B"]] # size of jumps
  C <- pars[["C"]] # size of sinusoidal waves
  mu <- pars[["mu"]] # constant
  P <- pars[["P"]] # period of jumps/waves
  
  x <- floor(t/(P*pi))
  y <- floor(t/(2*P*pi))
  
  return(mu*t + pi*B*P*y + B*(x - 2*y)*(t - pi*P*x) - P*C*(1 - cos(t/P)))
            
}

mu <- function(t, pars){
  B <- pars[["B"]] # size of jumps
  C <- pars[["C"]] # size of sinusoidal waves
  mu <- pars[["mu"]] # constant
  P <- pars[["P"]] # period of jumps/waves
  
  x <- floor(t/(P*pi))
  y <- floor(t/(2*P*pi))
  
  if (x == 2*y){
    return(mu - C*sin(t/P))
  } else {
    return(mu + B - C*sin(t/P))
  }
  
}


library(stelfi)
source("Simulate.R")
dyn.load(TMB::dynlib("weibull_hawkes_inhomog"))
reps <- 10
runs <- 20
pars = expand.grid(mu = c(0.5), B = c(1), C = c(0.25), P = c(10), beta = c(0.01, 0.1),
                   abratio = c(0.25, 0.75), k = c(0.5, (2/3), 1.5, 2), times=c(400,1000,2500,6250))

output <- matrix(0, nrow = reps, ncol = 15)
output[,1] <- pars$mu[i]; output[,2] <- pars$B[i]
output[,3] <- pars$C[i]; output[,4] <- pars$P[i]
output[,5] <- pars$beta[i]; output[,6] <- pars$abratio[i]; output[,7] <- pars$k[i]
output[,15] <- pars$times[i]

mu_max <- pars$mu[i] + pars$B[i] + abs(pars$C[i])
mu_min <- pars$mu[i] - abs(pars$C[i])
baseline_pars <- list(B = pars$B[i], C = pars$C[i], mu = pars$mu[i], P = pars$P[i])

incr <- 0.1
ts <- seq(incr, pars$times[i], by = incr)
covs <- -1*sin(ts/pars$P[i])
states <- rep(0, length = length(ts))

for(ind in 1:length(ts)){
  x <- floor(ts[ind]/(pars$P[i]*pi))
  y <- floor(ts[ind]/(2*pars$P[i]*pi))
  
  if (x > (2*y)){
    states[ind] = 1
  }
}

for(rep in 1:reps){
  times = simWIHP(Mu, mu_max, mu_min, pars$abratio[i]*pars$beta[i], pars$beta[i], 
                  baseline_pars, pars$k[i], pars$times[i], 1.e-9)$times
  NLLmin <- Inf
  slots <- floor(times/incr)
  for(j in 1:runs){
    mu <- runif(1, 0.1, 10)
    B <- runif(1, 0.1, 10)
    C <- runif(1, 0.1, 0.9) * mu
    beta <- runif(1, 0.1, 10)
    alpha <- runif(1, 0.1, 0.9) * beta
    k <- runif(1, 0.25, 4)
    out <- tryCatch(
      {
      obj <- TMB::MakeADFun(data = list(times = times, marks = rep(1, length(times)),
                                        slots = slots, covs = covs, states = states,
                                        interval = incr),
                        parameters = list(baseline = c(log(mu), log(B), qlogis(C/mu)),
                                            logit_abratio = qlogis(alpha/beta),
                                            log_beta = log(beta),
                                            log_k = log(k)),
                        hessian = TRUE, DLL = "weibull_hawkes_inhomog", silent = TRUE)
    
      opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = 0))
      obj$objective <- opt$value
      if (obj$objective < NLLmin && obj$objective != 0){
        z <- stelfi::get_coefs(obj)
        NLLmin <- obj$objective
      }
      
      output[rep, 8:10] <- z[1:3,1]; output[rep, 11] <- pars$P[i]
      output[rep, 13] <- z[4,1]/z[5,1]; output[rep, 12] <- z[5,1]
      output[rep, 14] <- z[6,1]
      }, error = function(cond){
        print(cond)
      }
    )
  }
}

write.csv(output, paste("Output/WIHPsims_",i,"_.csv",sep=""))
