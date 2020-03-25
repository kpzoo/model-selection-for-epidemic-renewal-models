######################################################################
## Negative-binomial posterior and APE score 
# From: Parag, KV, and Donnelly, CA. (2019) “Optimising Renewal Models for 
# Real-Time Epidemic Prediction and Estimation” BioRxiv: 835181.
######################################################################

# Assumptions
# - uses Poisson renewal equation (as in EpiEstim)
# - gamma prior required on R => negative binomial predictions

# Inputs - window length (k), SI distribution (sidistr), total infectiousness (Lday),
# incidence (Iday), max time (tPts), confidence 100*c(a, 1-a)%, gamma shape-scale (Rprior),
# start time for all windows (trunc)

# Output - APE score (ape), prob of next data (prob), R(t) estimate mean and confidence (Rhat, Rhatci), 
# I(t+1) prediction mean and confidence (Inexhat, Inexci)

apePredPost <- function(k, sidistr, Lday, Iday, Rprior, tmax, a, trunc){
  
  # Time over which estimation done
  idPts = seq(trunc, tmax-1); nPts = length(idPts)
  # Grouped incidence and infectiousness
  B = rep(0, nPts); A = B;
  
  # Offset Lday so for next point i.e. t+1
  Lday = Lday[2:length(Lday)] 
  
  # At each time before compute historical R
  for(i in 1:nPts){
    # Max window time, truncate if negative going backwards
    idmax = max(idPts[i] - k + 1, 1)
    # Look-back window of k (or less) indices
    idback = seq(idPts[i], idmax, -1) 
    # Relevant incidence sum (B) and total infectiousness sum (A)
    B[i] = sum(Iday[idback]); A[i] = sum(Lday[idback]);
  }

  # Shape-scale gamma R(t) parameters
  alpha = Rprior[1] + B
  beta = 1./(1/Rprior[2] + A)
  # Posterior predictive negative binomial parameter
  pr = Lday[idPts]*beta; pr = pr/(1 + pr)
  
  # Posterior mean estimate of R(t)
  Rhat = alpha*beta
  # Confidence interval (95%) on R(t) estimate
  Rhatci = matrix(-1, 2, length(Rhat))
  Rhatci[1,] = qgamma(a, shape = alpha, scale = beta)
  Rhatci[2,] = qgamma(1-a, shape = alpha, scale = beta)
  
  # Posterior mean prediction of I(t+1)
  Inexhat = Lday[idPts]*Rhat # check Lam[i] is for i+1 prediction <---
  #Inexhat = qnbinom(0.5, size = alpha, prob = 1-pr) 
  
  # Confidence interval (95%) on I(t+1) projection
  Inexci = matrix(-1, 2, length(Inexhat))
  Inexci[1,] = qnbinom(a, size = alpha, prob = 1-pr) 
  Inexci[2,] = qnbinom(1-a, size = alpha, prob = 1-pr) 
  
  # Probability of next incidence value at t+1s
  prob = dnbinom(Iday[seq(trunc+1, tmax)], size = alpha, prob = pr)
  # APE score at this k, assume log 1/0 = 0 <----- check
  ape = -sum(log(prob[prob != 0]))
  # PMSE score
  pmse = sum((Inexhat - Iday[seq(trunc+1, tmax)])^2)
  
  # Flag zero probabilities
  if(any(prob == 0)){
   # print(paste0(c('Zero prediction probs at:', k), collapse = ' '))
  }
  
  # Main outputs including parameters
  apePredPost = list(ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci, alpha, beta, pr)
}