# Fitting the Model expressed in Equations 4 and 8
# Fitted with NIMBLE in a Bayesian framework


library(nimble)

dWHPF <- nimbleFunction(
  run = function(x = double(2), # event times, event slots
                 N = integer(0),
                 n_slots = integer(0),
                 interval = double(0),
                 drs = double(2), # depths, rates, states
                 baseline = double(1),
                 ABratio = double(0),
                 beta = double(0),
                 k = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    times <- x[,1]
    slots <- x[,2]
    
    depths <- drs[,1]
    rates <- drs[,2]
    states <- drs[,3]
    
    alpha = ABratio * beta
    
    
    baselines <- numeric(n_slots)
    for(j in 1:n_slots){
      if (states[j] > 0){
        baselines[j] = exp(baseline[1] + baseline[2] + (baseline[3] * depths[j]) + (baseline[4] * rates[j]))
      } else {
        baselines[j] = exp(baseline[1] + (baseline[3] * depths[j]) + (baseline[4] * rates[j]))
      }
    }
    
    baselines_sums = numeric(n_slots)
    baselines_sums[1] = baselines[1]
    for(j in 2:n_slots){
      baselines_sums[j] = baselines_sums[(j-1)] + baselines[j]
    }
    
    A <- numeric(N)
    for(i in 2:N){
      # Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
      A[i] = exp(-beta * (times[i] - times[(i - 1)])) * (1.0 + A[(i - 1)])
    }
    
    compensators <- numeric(N)
    lambdas <- numeric(N)
    baseline_integral <- 0
    
    for(i in 1:N){
      slot = slots[i]
      lambdas[i] = baselines[slot] + alpha * A[i] # lambda at time i
      if (slot > 1){ # integral of baseline from 0 to time i
        baseline_integral = baselines_sums[(slot-1)] * interval
        baseline_integral = baseline_integral + (times[i] - interval*slot) * baselines[slot]
      } else {
        baseline_integral = times[i] * baselines[1]
      }
      compensators[i] = baseline_integral + (alpha/beta) * (i - 1 - A[i])
    }
    
    diff_compensators <- numeric(N)
    diff_compensators[1] = compensators[1]
    for(i in 2:N){
      diff_compensators[i] = compensators[i] - compensators[(i-1)]
    }
    # likelihood is density of standardized IATs (from a Weibull distribution)
    # multiplied by rates at each event
    
    nll = 0
    lgk = lgamma(1+(1/k)) # Avoid having to recompute
    gk = exp(lgk)
    
    nll = nll - sum(log(lambdas))
    nll = nll - N * log(k)
    nll = nll - N * k * lgk
    nll = nll - (k - 1) * sum(log(diff_compensators))
    for(i in 1:N){
      nll = nll + (diff_compensators[i] * gk) ^ k
    }
    
    if(log) return(-1*nll)
    else return(exp(-1*nll))
  })

assign('dWHPF', dWHPF, .GlobalEnv)

mcWHPF <- nimbleCode({
  # Priors
  ABratio ~ dunif(0,(1-1e-20))
  beta ~ dunif(1e-20,10)
  k ~ dunif(0.01,100)
  baseline[1] ~ dunif(-30,0)
  baseline[2] ~ dunif(0,100)
  baseline[3] ~ dunif(-2,2)
  baseline[4] ~ dunif(-100,100)
  
  
  # Likelihood
  events[,] ~ dWHPF(N = N, n_slots = n_slots, interval = interval, drs = drs[,], 
                    baseline = baseline[], ABratio = ABratio, beta = beta,
                    k = k)
  
  
})

## START OF SCRIPT
load("cleaned_whale_data.Rdata")
wh <- 53

times <- cleaned_data[[wh]]$times
depths <- cleaned_data[[wh]]$depths
rates <- cleaned_data[[wh]]$rates
incr <- cleaned_data[[wh]]$incr

drs <- matrix(0, nrow = length(depths), ncol = 3)
drs[,1] <- depths
drs[,2] <- rates
drs[,3] <- depths > 0.2


constants <- list(N = length(times), n_slots = length(depths), interval = incr, drs = drs)


initsList <- list(k = 1, ABratio = 0.5, beta = 1, baseline = c(-5,0,0,0))

events <- matrix(0, nrow = length(times), ncol = 2)
events[, 1] <- times
events[, 2] <- ceiling(times/incr)
data <- list(events = events)


# Run MCMC
Hmodel <- nimbleModel(mcWHPF, constants = constants, 
                      data = data, inits = initsList)

mcmc.out <- nimbleMCMC(model = Hmodel,
                       niter = 50000, nchains = 1, nburnin = 10000, thin =4, 
                       summary = TRUE, WAIC = TRUE)

# Save Results
write.csv(mcmc.out$samples, paste("Output/NIMBLE_Samples",wh,".csv", sep=""))
write.csv(mcmc.out$summary, paste("Output/NIMBLE_Coefs",wh,".csv", sep=""))
