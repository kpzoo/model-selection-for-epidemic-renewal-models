% APE for prediction and model selection in renewal models
function [apeIDs, pmseIDs, tday, Rhat, Rhatpmse, Ihat, Ihatpmse, ape, pmse, apeMod, pmseMod,...
    Rtrue, Iday, logPred, msePred] = apeForecastSwitchFn(ts, Rs, tday0, ks, nks, scenNo)

% Assumptions and notes
% - runs APE for single step switch (swtype 0 is downward)
% - includes APE across time as well (cumulatively)
% - uses NB posterior predictions and Bayesian APE
% - APE compares sequential predictions with true values
% - simulates a single epidemic, predicts I(t), estimates R(t)

% Simulate epidemic scenarios and truncate
Iwarn = 1; % ensure no warnings
while Iwarn 
    [Iday, Lam, Rtrue, tday, Iwarn] = epiSimSwitch(length(tday0), 1, ts, Rs, scenNo);
    % Ensure epidemic did not die out
    if length(find(Iday == 0)) > 50
        Iwarn = 1;
    end
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period based on removed 0s
nday = length(tday);

%% APE model selection and prediction

% Posterior I predictions and R estimates
pred = cell(1, nks); R = pred; 
% APE metric and PMSE
prob = pred; ape = zeros(1, nks); pmse = ape;

for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, prob{i}, R{i}, ~] = getNegBinSwQuick(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
end

% Predictive (successive) prob  
probS = cell2mat(prob');
logPred = log(1./probS);
logPred = cumsum(logPred');
% Successive APE with time (local)
[~, apeIDs] = min(logPred, [], 2); 
apeIDs = apeIDs';

% Mean prediction error with time
predS = cell2mat(pred');
msePred = (predS - Iday(2:end)).^2;
% Sum of squares then mean
msePred = cumsum(msePred');
for i = 1:nday-1
    msePred(i, :) = msePred(i, :)/i;
end
[~, pmseIDs] = min(msePred, [], 2);
pmseIDs = pmseIDs';

% Best models according to metrics
[~, apeMod] = min(ape);
[~, pmseMod] = min(pmse);

% Best mean R estimate and mean I prediction
Rhat = R{apeMod}; Rhatpmse = R{pmseMod};
Ihat = pred{apeMod}; Ihatpmse = pred{pmseMod};

