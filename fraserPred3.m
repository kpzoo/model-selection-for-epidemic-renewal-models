% APE for prediction and model selection in renewal models
clearvars;
clc; close all;

% Assumptions and notes
% - Directly compares the Bayesian and MLE results
% - Posterior predictions included (negative binomial Bayesian)
% - APE measures with actual value and its support in predictive distrib
% - simulates a single epidemic under Fraser-Cori renewal model

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Time code
tic;
% Save data
saveTrue = 0;

% Define a scenario
scenNo = 2;
% Time for epidemic observation (days)
tday = 1:401; nday = length(tday);
% Simulate epidemic scenarios
[Iday, Lam, Rtrue, Ifil, tday, Iwarn] = epiSimRenew(scenNo, tday, nday, 0, 1);
lenRt = length(Rtrue);

%% APE model selection and prediction

% Parameters for model selection
ks = unique(round(linspace(2, 200, 100))); 
nks = length(ks);
disp(['k varies from ' num2str(ks(1)) ' to ' num2str(ks(end))]);

% Predictions at each k (MLE based)
pred = cell(1, nks); predInt = pred; prob = pred; 
% Bayesian predictions
predB = pred; predIntB = predInt; probB = prob; CV = pred; predB2 = pred;
% Metrics for model selection
ape = zeros(1, nks); pmse = ape; apeB = ape; pmseB = ape;

for i = 1:nks
    % One step ahead MLE predictions at given k
    [pred{i}, predInt{i}, prob{i}, ~] = getAPE(ks(i), lenRt, Iday, Lam);
    % One step ahead Bayesian posterior prediction for k
    [predB{i}, predIntB{i}, probB{i}, ~, ~] = getNegBin(ks(i), lenRt, Iday, Lam);
    
    % APE and predictive MSE for MLE
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
    % APE and predictive MSE for Bayesian
    apeB(i) = -sum(log(probB{i}));
    pmseB(i) = mean((Iday(2:end) - predB{i}).^2);
    
    disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
end

% Best models according to metrics
[apeMin, apeMod] = min(ape);
[pmseMin, pmseMod] = min(pmse);
[apeMinB, apeModB] = min(apeB);
[pmseMinB, pmseModB] = min(pmseB);

% Best ks and nGrps
modID = [apeMod pmseMod];
kbest = ks(modID); 
disp(['MLE k: [ape pmse] = [' num2str(kbest) ']' ]);
modIDB = [apeModB pmseModB];
kbestB = ks(modIDB); 
disp(['Bayes k: [ape pmse] = [' num2str(kbestB) ']' ]);

%% Visualisation and post processing

% Examine predictions against true incidence MLE
figure;
subplot(2, 2, 1);
plot(tday(2:end), Iday(2:end), 'c', 'linewidth', 2);
hold on;
plot(tday(2:end), pred{apeMod}, 'color', grey2, 'linewidth', 2);
plot(tday(2:end), predInt{apeMod}(1, :), 'color', grey1, 'linewidth', 2);
plot(tday(2:end), predInt{apeMod}(2, :), 'color', grey1, 'linewidth', 2);
h = gca; h.XLim = [tday(2) tday(end)];
if max(Iday) > 2000
    h.YScale = 'log';
end
legend(['k = ' num2str(kbest(1))], 'location', 'best');
hold off; grid off; box off;
xlabel('time (days)');
ylabel('APE MLE');
subplot(2, 2, 2);
plot(tday(2:end), Iday(2:end), 'c', 'linewidth', 2);
hold on;
plot(tday(2:end), pred{pmseMod}, 'color', grey2, 'linewidth', 2);
plot(tday(2:end), predInt{pmseMod}(1, :), 'color', grey1, 'linewidth', 2);
plot(tday(2:end), predInt{pmseMod}(2, :), 'color', grey1, 'linewidth', 2);
h = gca; h.XLim = [tday(2) tday(end)];
if max(Iday) > 2000
    h.YScale = 'log';
end
legend(['k = ' num2str(kbest(2))], 'location', 'best');
hold off; grid off; box off;
xlabel('time (days)');
ylabel('PMSE MLE');

% Examine predictions against true incidence Bayes
subplot(2, 2, 3);
plot(tday(2:end), Iday(2:end), 'c', 'linewidth', 2);
hold on;
plot(tday(2:end), predB{apeModB}, 'color', grey2, 'linewidth', 2);
plot(tday(2:end), predIntB{apeModB}(1, :), 'color', grey1, 'linewidth', 2);
plot(tday(2:end), predIntB{apeModB}(2, :), 'color', grey1, 'linewidth', 2);
h = gca; h.XLim = [tday(2) tday(end)];
if max(Iday) > 2000
    h.YScale = 'log';
end
legend(['k = ' num2str(kbestB(1))], 'location', 'best');
hold off; grid off; box off;
xlabel('time (days)');
ylabel('APE Bayes');
subplot(2, 2, 4);
plot(tday(2:end), Iday(2:end), 'c', 'linewidth', 2);
hold on;
plot(tday(2:end), predB{pmseModB}, 'color', grey2, 'linewidth', 2);
plot(tday(2:end), predIntB{pmseModB}(1, :), 'color', grey1, 'linewidth', 2);
plot(tday(2:end), predIntB{pmseModB}(2, :), 'color', grey1, 'linewidth', 2);
h = gca; h.XLim = [tday(2) tday(end)];
if max(Iday) > 2000
    h.YScale = 'log';
end
legend(['k = ' num2str(kbestB(2))], 'location', 'best');
hold off; grid off; box off;
xlabel('time (days)');
ylabel('PMSE Bayes');

% Model selection with APE and PMSE
figure;
yyaxis left
plot(ks, ape, 'c--', 'linewidth', 2);
hold on;
plot(ks, apeB, 'b--', 'linewidth', 2);
hold off;
ylabel('APE');
yyaxis right
plot(ks, pmse, 'g', 'linewidth', 2);
hold on;
plot(ks, pmseB, 'r--', 'linewidth', 2);
hold off;
ylabel('PMSE');
grid off; box off;
xlabel('$k$ (days)');
legend('APE MLE', 'APE Bayes', 'PMSE MLE', 'PMSE Bayes', 'location', 'best');


%% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);


