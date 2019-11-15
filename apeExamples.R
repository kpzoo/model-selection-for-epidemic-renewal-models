# Extract incidence data from EpiEstim package

# Assumptions and modifications
# - incidence data for APE calculation in matlab
# - obtain Cori R estimates and output

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
trajName = 'EpiEstim test'
dir.create(file.path(this.dir, trajName))
pathf = paste(c(this.dir, '/', trajName, '/'), collapse = "")

# Load data on pandemic flu in Baltimore in 1918
data("Flu1918");
Iflu = Flu1918$Incidence
genflu = Flu1918$SI.Distr
# Total infectiousness
Lflu = OverallInfectivity(Iflu, genflu)

# Estimates on a 7 day window - considered best in Cori 2013 for this data
estflu = EstimateR(Iflu, T.Start=2:86, T.End=8:92, method="NonParametricSI", SI.Distr = genflu)
#estflu = EstimateR(Iflu, T.Start=2:91, T.End=3:92, method="NonParametricSI", SI.Distr = genflu)
tflu = estflu$R$T.End # end of window
Rflu = estflu$R$`Mean(R)`
RfluCI = matrix(NA, nrow = 2, ncol = length(Rflu))
RfluCI[1,] = estflu$R$`Quantile.0.025(R)`
RfluCI[2,] = estflu$R$`Quantile.0.975(R)`

# Load data on SARS in Hong Kong in 2003
data("SARS2003")
Isars = SARS2003$Incidence
gensars = SARS2003$SI.Distr
# Total infectiousness
Lsars = OverallInfectivity(Isars, gensars)

# Estimates on a 7 day window - considered best in Cori 2013 for this data
estsars = EstimateR(Isars, T.Start=14:101, T.End=20:107, method="NonParametricSI", SI.Distr = gensars)
#estsars = EstimateR(Isars, T.Start=14:106, T.End=15:107, method="NonParametricSI", SI.Distr = gensars)
tsars = estsars$R$T.End # end of window
Rsars = estsars$R$`Mean(R)`
RsarsCI = matrix(NA, nrow = 2, ncol = length(Rsars))
RsarsCI[1,] = estsars$R$`Quantile.0.025(R)`
RsarsCI[2,] = estsars$R$`Quantile.0.975(R)`

# Output all data
tableWrite(Iflu, 'Iflu.csv', pathf); tableWrite(Isars, 'Isars.csv', pathf)
tableWrite(Lflu, 'Lflu.csv', pathf); tableWrite(Lsars, 'Lsars.csv', pathf)
tableWrite(genflu, 'genflu.csv', pathf); tableWrite(gensars, 'gensars.csv', pathf)
tableWrite(tflu, 'tflu.csv', pathf); tableWrite(tsars, 'tsars.csv', pathf)
tableWrite(Rflu, 'Rflu.csv', pathf); tableWrite(Rsars, 'Rsars.csv', pathf)
tableWrite(RfluCI, 'RCIflu.csv', pathf); tableWrite(RsarsCI, 'RCIsars.csv', pathf)

# Estimate R at each day over 7 day window finishing on that day
quartz()
EstimateR(Flu1918$Incidence, T.Start=2:86, T.End=8:92, method="NonParametricSI",
          SI.Distr=Flu1918$SI.Distr, plot=TRUE, leg.pos=xy.coords(60,2.5))
dev.copy2eps(file="fluCori.eps")
quartz()
EstimateR(SARS2003$Incidence, T.Start=14:101, T.End=20:107, method="NonParametricSI",
          SI.Distr=SARS2003$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,7))
dev.copy2eps(file="sarsCori.eps")
