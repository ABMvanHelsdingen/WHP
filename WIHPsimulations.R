# simulate WHP with inhomogenous baseline rate
# Baseline rate is mu + B * s(t) + C * sin(t)
# where B, C are non-negative and mu > C
# and s(t) is 1 when pi<t<2*pi, 3*pi<t<4*pi, 5*pi<t<6*pi etc

# following block allows the script to be run in parallel on NeSI
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
set.seed(4321 + 10*i)

# Integral of baseline

Mu <- function(t, pars){
  B <- pars[["B"]]
  C <- pars[["C"]]
  mu <- pars[["mu"]]
  
  x <- floor(t/pi)
  y <- floor(t/(2*pi))
  
  return(mu*t + pi*B*y + B*(x - 2*y)*(t - pi*x) + C*(1 - cos(t)))

  
            
}


library(stelfi)
source("Ogata.R")
dyn.load(TMB::dynlib("weibull_hawkes_inhomog"))
reps <- 100
runs <- 20
pars = expand.grid(mu = c(1), B = c(1), C = c(0.5), beta = c(0.01, 0.1),
                   abratio = c(0.25, 0.75), k = c(0.5, (2/3), 1.5, 2), times=c(400,1000,2500,6250))

output <- matrix(0, nrow = reps, ncol = 13)
output[,1] <- pars$mu[i]; output[,2] <- pars$A[i]
output[,3] <- pars$B[i]; output[,4] <- pars$beta[i]
output[,5] <- pars$abratio[i]; output[,6] <- pars$k[i]
output[,13] <- pars$times[i]

mu_max <- pars$mu[i] + pars$B[i] + pars$C[i]
mu_min <- pars$mu[i] - pars$C[i]
baseline_pars <- list(B = pars$B[i], C = pars$C[i], mu = pars$mu[i])

incr <- 0.1
ts <- seq(incr, pars$times[i], by = incr)
covs <- sin(ts)
states <- rep(0, length = length(ts))

for(ind in 1:length(ts)){
  x <- floor(ts[ind]/pi)
  y <- floor(ts[ind]/(2*pi))
  
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
      
      output[rep, 5] <- z[1,1]; output[rep, 6] <- z[3,1]
      output[rep, 7] <- z[2,1]/z[3,1]; output[rep, 8] <- z[4,1]
      }, error = function(cond){
        print(cond)
      }
    )
  }
}

write.csv(output, paste("Output/WIHPsims_",i,"_.csv",sep=""))

