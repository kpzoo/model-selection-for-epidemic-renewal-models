######################################################################
## Plot APE estimate and prediction
# From: Parag, KV, and Donnelly, CA. (2019) “Optimising Renewal Models for 
# Real-Time Epidemic Prediction and Estimation” BioRxiv: 835181.
######################################################################

# Assumptions
# - uses model output (Rmod) from apeEstim.R as [[ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci]]

# Inputs -  apeEstim.R list i.e. Rmod[[2]] for APE or Rmod[[4]] for PMSE, plot name string,
# best window length so Rmod[[1]] for APE and Rmod[[3]] for PMSE, true incidence Iplt
# Output - .eps plots of best R estimate and one-step-ahead incidence predictions

plotAPEWindow <- function(Rmod, plotname, kbest, Iplt){
  # Extract best I(t+1) and R(t) estimates/predictions
  Rhat = Rmod[[4]]; Rhatci = Rmod[[5]]
  Inexhat = Rmod[[6]]; Inexhatci = Rmod[[7]]
  
  # Check lengths
  if (length(Rhat) != length(Inexhat)){
    print('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    tset = 1:length(Rhat)
    
    # Plot as composite
    #setEPS()
    #postscript(paste0(c(plotname, '.eps'), collapse = ''))
    
    quartz()
    par(mfrow=c(2,1))
    # Reprod. num estimates and confidence interval
    plot(tset, Rhat, type = 'l', bty = 'l', lwd = 2, col='red',
         xlab = paste0("time (k = ", kbest, ")"), ylab = 'Rhat')
    lines(tset, Rhatci[1,], col = 'red', type = "l", lwd = 1)
    lines(tset, Rhatci[2,], col = 'red', type = "l", lwd = 1)
    # plot(tset, Rhat, pch = 19, bty = 'l', lwd = 2, col='red',
    #      xlab = paste0("time (k = ", kbest, ")"), ylab = 'Rhat')
    # arrows(x0 = tset, y0 = Rhatci[1,], x1 = tset, y1 = Rhatci[2,], code=3, angle=90, length=0.1, col = 'lightgray')
    
    # Incidence predictions and confidence interval
    plot(tset, Inexhat, type = 'l', bty = 'l', lwd = 2, col='red',
         xlab = paste0("time (k = ", kbest, ")"), ylab = 'Ihat')
    lines(tset, Inexhatci[1,], col = 'red', type = "l", lwd = 1)
    lines(tset, Inexhatci[2,], col = 'red', type = "l", lwd = 1)
    points(tset, Iplt, pch = 19, col = 'gray')
    
    dev.copy2eps(file=paste0(c(plotname, '.eps'), collapse = ''))
  }
  
  
}