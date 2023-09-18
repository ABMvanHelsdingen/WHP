# Fitting the Model expressed in Equations 9 and 10
# Fitted with TMB in a frequentist framework

load("cleaned_whale_data.Rdata")
wh <- 53
dyn.load(TMB::dynlib("whale_cues"))
library(stelfi)


n_starts <- 200
cbfNLLs <- matrix(0, nrow = 1, ncol = n_starts)
cbfCoefs <- matrix(0, nrow = n_starts, ncol = 21)

times <- cleaned_data[[wh]]$times
depths <- cleaned_data[[wh]]$depths
rates <- cleaned_data[[wh]]$rates
rincr <- 1/cleaned_data[[wh]]$incr # reciprocal of the increment
  
# assign states
states <- numeric(length(depths))
for(j in 1:length(depths)){
  if (depths[j] > 0.2 && rates[j] >= 0){
    states[j] <- 1 # descent
  } else if (depths[j] > 0.2 && rates[j] < 0){
    states[j] <- 2 # ascent
  } # if depth < 20m, state = 0 (surface)
}
  
  

# calculate slot number of each event
slot_numbers <- floor(times * rincr)
NLLmin <- Inf
for(j in 1:n_starts){
  out <- tryCatch(
  {
    rbeta <- rexp(1,1)
    pars <- list(alpha = runif(1,0.01,0.99) * rbeta ,beta = rbeta, k = runif(1, 0.8, 5))
    background_pars <- c(runif(1,-12,-6), runif(1,-2,0), runif(1,-1,1), runif(1,-7,-0.5))
    cbfCoefs[j, 1:7] <- c(background_pars[1], background_pars[2], background_pars[3], background_pars[4],
                          pars[[1]], pars[[2]], pars[[3]])
    obj <- TMB::MakeADFun(data = list(times = times, marks = rep(1, length(times)), interval = 1/ rincr,
                                      slots = slot_numbers, states = states, depths = depths, rates = rates),
                          parameters = list(baseline = background_pars,
                                            logit_abratio = qlogis(pars[["alpha"]]/pars[["beta"]]),
                                            log_beta = log(pars[["beta"]]),
                                            log_k = log(pars[["k"]])),
                          hessian = TRUE, DLL = "whale_cues", silent = TRUE)
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = 0))
    obj$objective <- opt$value

    cbfNLLs[1, j] <- obj$objective
    coefs <- stelfi::get_coefs(obj)
    cbfCoefs[j, 8:14] <- coefs[,1] # MLEs
    cbfCoefs[j, 15:21] <- coefs[,2] # MLE standard errors

    if ((obj$objective < NLLmin) && (obj$objective != 0)) {
      NLLmin <- obj$objective
    }
  }, error = function(cond){
    print(paste("Fitting error", wh, j))
    print(cond)
  }
  )
}

# Save Results
write.csv(cbfCoefs,paste("Output/FreqCoefs_",wh, ".csv", sep=""))
write.csv(cbfNLLs,paste("Output/FreqNLLs_",wh, ".csv", sep=""))