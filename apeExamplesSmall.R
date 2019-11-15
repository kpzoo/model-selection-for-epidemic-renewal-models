# Extract incidence data from EpiEstim package

# Assumptions and modifications
# - incidence data for APE calculation in matlab
# - obtain Cori R estimates and output
# - only small epidemics considered 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main packages
library("EpiEstim")

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Function to write simple csv files to correct path
tableWrite <- function(val, name, pathname) {
  # Add path to name
  str0 <- paste(c(pathname, name), collapse = "")
  # Write table
  write.table(val, str0, row.names=FALSE, col.names=FALSE, sep=",")
}

# Create folder for traj specific results
trajName = 'EpiEstim data'
dir.create(file.path(this.dir, trajName))
pathf = paste(c(this.dir, '/', trajName, '/'), collapse = "")

# Load data on measles in Hallegoch in 1861
data("Measles1861")
Imeas = Measles1861$Incidence
genmeas = Measles1861$SI.Distr
# Total infectiousness
Lmeas = OverallInfectivity(Imeas, genmeas)

# Estimates on a 7 day window - considered best in Cori 2013 for this data
estmeas = EstimateR(Imeas,  T.Start=17:42, T.End=23:48, method="NonParametricSI", SI.Distr = genmeas)
tmeas = estmeas$R$T.End # end of window
Rmeas = estmeas$R$`Mean(R)`
RmeasCI = matrix(NA, nrow = 2, ncol = length(Rmeas))
RmeasCI[1,] = estmeas$R$`Quantile.0.025(R)`
RmeasCI[2,] = estmeas$R$`Quantile.0.975(R)`

# Load data on on smallpox in Kosovo in 1972
data("Smallpox1972")
Ipox = Smallpox1972$Incidence
genpox = Smallpox1972$SI.Distr
# Total infectiousness
Lpox = OverallInfectivity(Ipox, genpox)

# Estimates on a 7 day window - considered best in Cori 2013 for this data
estpox = EstimateR(Ipox, T.Start=27:51, T.End=33:57, method="NonParametricSI", SI.Distr = genpox)
tpox = estpox$R$T.End # end of window
Rpox = estpox$R$`Mean(R)`
RpoxCI = matrix(NA, nrow = 2, ncol = length(Rpox))
RpoxCI[1,] = estpox$R$`Quantile.0.025(R)`
RpoxCI[2,] = estpox$R$`Quantile.0.975(R)`

# Output all data
tableWrite(Imeas, 'Imeas.csv', pathf); tableWrite(Ipox, 'Ipox.csv', pathf)
tableWrite(Lmeas, 'Lmeas.csv', pathf); tableWrite(Lpox, 'Lpox.csv', pathf)
tableWrite(genmeas, 'genmeas.csv', pathf); tableWrite(genpox, 'genpox.csv', pathf)
tableWrite(tmeas, 'tmeas.csv', pathf); tableWrite(tpox, 'tpox.csv', pathf)
tableWrite(Rmeas, 'Rmeas.csv', pathf); tableWrite(Rpox, 'Rpox.csv', pathf)
tableWrite(RmeasCI, 'RCImeas.csv', pathf); tableWrite(RpoxCI, 'RCIpox.csv', pathf)

# Estimate R at each day over 7 day window finishing on that day
quartz()
EstimateR(Measles1861$Incidence, T.Start=17:42, T.End=23:48, method="NonParametricSI",
          SI.Distr=Measles1861$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,7))
dev.copy2eps(file="measCori.eps")
quartz()
EstimateR(Smallpox1972$Incidence, T.Start=27:51, T.End=33:57, method="NonParametricSI",
          SI.Distr=Smallpox1972$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,15))
dev.copy2eps(file="poxCori.eps")
