######################################################################
## Demonstrate APE with EpiEstim package
######################################################################

# Assumptions and modifications
# - compute EpiEstim for influenza and SARS data
# - get APE-regularised estimates for comparison

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main packages
library("EpiEstim")
library("caTools")

# Key functions
source('apeEstim.R')
source('apePredPost.R')
source('apeSpecific.R')
source('plotAPEWindow.R')

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Boolean for plotting
wantPlot = 0

# Load data on pandemic flu in Baltimore in 1918
data("Flu1918");
Iflu = Flu1918$incidence
genflu = Flu1918$si_distr
# Total infectiousness
Lflu = overall_infectivity(Iflu, genflu)

# Load data on SARS in Hong Kong in 2003
data("SARS2003")
Isars = SARS2003$incidence
gensars = SARS2003$si_distr
# Total infectiousness
Lsars = overall_infectivity(Isars, gensars)

######################################################################
## Conventional EpiEstim estimates
######################################################################

# Start and end times for weekly window in flu
t_start = seq(2, length(Iflu)-6)
t_end = t_start + 6
# Estimates on a 7 day window - considered best in Cori 2013 for this data
estflu = estimate_R(Iflu, method ="non_parametric_si", config = make_config(list(
  si_distr = genflu, t_start = t_start, t_end = t_end)))
# Extract outputs
tflu = estflu$R$t_end # end of window
Rflu = estflu$R$`Mean(R)`
RfluCI = matrix(NA, nrow = 2, ncol = length(Rflu))
RfluCI[1,] = estflu$R$`Quantile.0.025(R)`
RfluCI[2,] = estflu$R$`Quantile.0.975(R)`

# Start and end times for weekly window in SARS
t_start = seq(14, length(Isars)-6)
t_end = t_start + 6
# Estimates on a 7 day window - considered best in Cori 2013 for this data
estsars = estimate_R(Isars, method ="non_parametric_si", config = make_config(list(
  si_distr = gensars, t_start = t_start, t_end = t_end)))
# Extract outputs
tsars = estsars$R$t_end # end of window
Rsars = estsars$R$`Mean(R)`
RsarsCI = matrix(NA, nrow = 2, ncol = length(Rsars))
RsarsCI[1,] = estsars$R$`Quantile.0.025(R)`
RsarsCI[2,] = estsars$R$`Quantile.0.975(R)`

# Estimate R at each day over 7 day window finishing on that day
if(wantPlot){
  quartz()
  plot(estflu)
  dev.copy2eps(file="fluCori.eps")
  quartz()
  plot(estsars)
  dev.copy2eps(file="sarsCori.eps")
}

######################################################################
## APE and PMSE solution 
######################################################################
# Outputs take list form [[kbest1, modBest1, kbest2, modBest2]]
# Each modBest is also a list of form [[ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci, alpha, beta, pr]]

# Priors and settings
Rprior = c(1, 5); a = 0.025; trunctime = 2

# Flu results
Rmodflu = apeEstim(Iflu, genflu, Lflu, Rprior, a, trunctime, 'flu')
# Best estimates and prediction 
plotAPEWindow(Rmodflu[[2]], 'flu', Rmodflu[[1]], Iplt = Iflu[seq(trunctime+1,length(Iflu))])
# Specific 7-day window
Rmodflu7 = apeSpecific(Iflu, genflu, Lflu, Rprior, a, trunctime, 'flu7', 7)
plotAPEWindow(Rmodflu7[[2]], 'flu7', Rmodflu7[[1]], Iplt = Iflu[seq(trunctime+1,length(Iflu))])

# Sars results
Rmodsars = apeEstim(Isars, gensars, Lsars, Rprior, a, trunctime, 'sars')
# Best estimates and prediction 
plotAPEWindow(Rmodsars[[2]], 'sars', Rmodsars[[1]], Iplt = Isars[seq(trunctime+1,length(Isars))])
# Specific 7-day window
Rmodsars7 = apeSpecific(Isars, gensars, Lsars, Rprior, a, trunctime, 'sars7', 7)
plotAPEWindow(Rmodsars7[[2]], 'sars7', Rmodsars7[[1]], Iplt = Isars[seq(trunctime+1,length(Isars))])