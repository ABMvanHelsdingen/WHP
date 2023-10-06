# simulate IWP
# following block allows the script to be run in parallel on NeSI
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
set.seed(4321 + 10*i)


library(stelfi)
source("Ogata.R")
dyn.load(TMB::dynlib("weibull_hawkes"))
reps <- 100
runs <- 20
pars = expand.grid(mu = c(1), beta = c(0.01, 0.1),
                   abratio = c(0.25, 0.75), k = c(0.5, (2/3), 1.5, 2), times=c(400,1000,2500,6250))

output <- matrix(0, nrow = reps, ncol = 9)
output[,1] <- pars$mu[i]; output[,2] <- pars$beta[i]
output[,3] <- pars$abratio[i]; output[,4] <- pars$k[i]
output[,9] <- pars$times[i]

for(rep in 1:reps){
  times = UWOgata(pars$mu[i], pars$abratio[i]*pars$beta[i], pars$beta[i], pars$k[i], pars$times[i])
  NLLmin <- Inf
  for(j in 1:runs){
    mu <- runif(1, 0.1, 10)
    beta <- runif(1, 0.1, 10)
    alpha <- runif(1, 0.1, 0.9) * beta
    k <- runif(1, 0.25, 4)
    out <- tryCatch(
      {
      obj <- TMB::MakeADFun(data = list(times = times, marks = rep(1, length(times))),
                        parameters = list(log_mu = log(mu),
                                            logit_abratio = qlogis(alpha/beta),
                                            log_beta = log(beta),
                                            log_k = log(k)),
                        hessian = TRUE, DLL = "weibull_hawkes", silent = TRUE)
    
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

write.csv(output, paste("Output/IWPsims_",i,"_.csv",sep=""))

