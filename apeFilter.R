# Extract incidence data from EpiEstim package

# Assumptions and modifications
# - updated to use new version of EpiEstim
# - also works on filtered incidence curves
# - incidence data for APE calculation in matlab
# - obtain Cori R estimates and output

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main packages
library("EpiEstim")
library("caTools")

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
trajName = 'EpiEstim filter'
dir.create(file.path(this.dir, trajName))
pathf = paste(c(this.dir, '/', trajName, '/'), collapse = "")

# Load data on pandemic flu in Baltimore in 1918
data("Flu1918");
Iflu = Flu1918$incidence
genflu = Flu1918$si_distr
# Total infectiousness
Lflu = overall_infectivity(Iflu, genflu)

# Start and end times for weekly window
t_start = seq(2, length(Iflu)-6)
t_end = t_start + 6
# Estimates on a 7 day window - considered best in Cori 2013 for this data
#estflu = EstimateR(Iflu, T.Start=2:86, T.End=8:92, method="NonParametricSI", SI.Distr = genflu)
#estflu = EstimateR(Iflu, T.Start=2:91, T.End=3:92, method="NonParametricSI", SI.Distr = genflu)
estflu = estimate_R(Iflu, method ="non_parametric_si", config = make_config(list(
  si_distr = genflu, t_start = t_start, t_end = t_end)))
tflu = estflu$R$t_end # end of window
Rflu = estflu$R$`Mean(R)`
RfluCI = matrix(NA, nrow = 2, ncol = length(Rflu))
RfluCI[1,] = estflu$R$`Quantile.0.025(R)`
RfluCI[2,] = estflu$R$`Quantile.0.975(R)`

# Load data on SARS in Hong Kong in 2003
data("SARS2003")
Isars = SARS2003$incidence
gensars = SARS2003$si_distr
# Total infectiousness
Lsars = overall_infectivity(Isars, gensars)

# Start and end times for weekly window
t_start = seq(14, length(Isars)-6)
t_end = t_start + 6

# Estimates on a 7 day window - considered best in Cori 2013 for this data
#estsars = EstimateR(Isars, T.Start=14:101, T.End=20:107, method="NonParametricSI", SI.Distr = gensars)
#estsars = EstimateR(Isars, T.Start=14:106, T.End=15:107, method="NonParametricSI", SI.Distr = gensars)
estsars = estimate_R(Isars, method ="non_parametric_si", config = make_config(list(
  si_distr = gensars, t_start = t_start, t_end = t_end)))
tsars = estsars$R$t_end # end of window
Rsars = estsars$R$`Mean(R)`
RsarsCI = matrix(NA, nrow = 2, ncol = length(Rsars))
RsarsCI[1,] = estsars$R$`Quantile.0.025(R)`
RsarsCI[2,] = estsars$R$`Quantile.0.975(R)`

# Compute moving averages over m points
m = 5
IfluFilt = runmean(Iflu, m)
IsarsFilt = runmean(Isars, m)
# Lam based on smoothed incidence curves
LfluFilt = overall_infectivity(IfluFilt, genflu)
LsarsFilt = overall_infectivity(IsarsFilt, gensars)
# Make integer
IfluFilt = round(IfluFilt); IsarsFilt = round(IsarsFilt)


# Compute median filters
mmed = 5
IfluFilt2 = runmed(Iflu, mmed); IfluFilt2 = as.numeric(IfluFilt2)
IsarsFilt2 = runmed(Isars, mmed); IsarsFilt2 = as.numeric(IsarsFilt2)
# Lam based on smoothed incidence curves
LfluFilt2 = overall_infectivity(IfluFilt2, genflu)
LsarsFilt2 = overall_infectivity(IsarsFilt2, gensars)
# Make integer
IfluFilt = round(IfluFilt2); IsarsFilt = round(IsarsFilt2)

# Output all data
tableWrite(Iflu, 'Iflu.csv', pathf); tableWrite(Isars, 'Isars.csv', pathf)
tableWrite(Lflu, 'Lflu.csv', pathf); tableWrite(Lsars, 'Lsars.csv', pathf)
tableWrite(genflu, 'genflu.csv', pathf); tableWrite(gensars, 'gensars.csv', pathf)
tableWrite(tflu, 'tflu.csv', pathf); tableWrite(tsars, 'tsars.csv', pathf)
tableWrite(Rflu, 'Rflu.csv', pathf); tableWrite(Rsars, 'Rsars.csv', pathf)
tableWrite(RfluCI, 'RCIflu.csv', pathf); tableWrite(RsarsCI, 'RCIsars.csv', pathf)
tableWrite(IfluFilt, 'IfluFilt.csv', pathf); tableWrite(IsarsFilt, 'IsarsFilt.csv', pathf)
tableWrite(LfluFilt, 'LfluFilt.csv', pathf); tableWrite(LsarsFilt, 'LsarsFilt.csv', pathf)
tableWrite(m, 'mfil.csv', pathf); tableWrite(mmed, 'mfil2.csv', pathf)
tableWrite(IfluFilt2, 'IfluFilt2.csv', pathf); tableWrite(IsarsFilt2, 'IsarsFilt2.csv', pathf)
tableWrite(LfluFilt2, 'LfluFilt2.csv', pathf); tableWrite(LsarsFilt2, 'LsarsFilt2.csv', pathf)

# Estimate R at each day over 7 day window finishing on that day
quartz()
plot(estflu)
dev.copy2eps(file="fluCori.eps")
quartz()
plot(estsars)
dev.copy2eps(file="sarsCori.eps")
