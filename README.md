# model-selection-for-epidemic-renewal-models
Uses accumulated prediction error (APE) to select among renewal model formulations 

Renewal models [1, 2] estimate the effective reproduction number, R(t), underyling observed outbreaks. They use incidence data (time series of infected cases) to infer a piecewise-constant approximation to R(t). The quality of this estimate is highly dependent on the size of a smoothing window (k) that is employed. This code presents a method for optimally selecting k in a manner that balances reliable R(t) estimation with short-term forecasts of incidence. This method is based on the accumulated prediction error (APE) idea from information theory [3].

The main code is in the m files denoted ape_ or batch_ and examines various epidemic examples using APE. The m files draw on data from the EpiEstim package in R [2] and are responsible for all figures and results in [4]. Figures and results are in the accompanying folders. An additional R implementation is also provided as APE_R.

[1] C. Fraser, D. Cummings, D. Klinkenberg, et al., “Influenza transmission in households during the 1918 pandemic,” Am. J. Epidemiol, vol. 174, no. 5, pp. 505–14, 2011.
[2] A. Cori, N. Ferguson, C. Fraser, et al., “A new framework and software to estimate time-varying reproduction numbers during epidemics,” Am. J. Epidemiol, vol. 178, no. 9, pp. 1505–12, 2013.
[3] J. Rissanen, “Order estimation by accumulated prediction errors,” J. Appl. Prob, vol. 23, pp. 55–61, 1986.
[4] K. Parag and C. Donnelly, "Optimising Renewal Models for Real-Time Epidemic Prediction and Estimation, " BioRxiv, 2019.
