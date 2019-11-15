% APE for renewal model prediction and selection with 2 steps ahead
clearvars; clc;
close all; tic;

% Assumptions and notes
% - computes 2 step prediction with posterior mean of I(t+1)
% - includes APE across time as well (cumulatively)
% - uses NB posterior predictions and Bayesian APE
% - APE compares sequential predictions with true values
% - simulates a single epidemic, predicts I(t), estimates R(t)

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving
saveFol = 'single data'; thisDir = cd;

% Define a scenario
scenNam = {'constant', 'cyclic', 'logistic', 'piecewise', 'boom-bust', 'bottle', '2-step', 'filtered'};
scenNo = 8; scenChoice = scenNam{scenNo};

% Inital time for epidemic observation (days)
tday = 1:201;
% Simulate epidemic scenarios and truncate
[Iday, Lam, Rtrue, tday, Iwarn] = epiSimAPE(scenNo, tday, length(tday), 1);
% Truncated observation period
nday = length(tday);

%% APE model selection and prediction

% Space of look-back windows
%ks = unique(round(linspace(2, floor(nday), 100)));
%ks = 1:nday;
ks = 1:ceil(nday/2);
nks = length(ks);
disp(['k varies from ' num2str(ks(1)) ' to ' num2str(ks(end))]);

% Posterior incidence predictions
pred = cell(1, nks); predInt = pred; pred2 = pred; 
predInt2 = pred; prob = pred; prob2 = pred;
% Posterior estimates of R over ks 
R = pred; RInt = pred; R2 = pred; RInt2 = pred;
% APE metric and PMSE
ape = zeros(1, nks); pmse = ape; ape2 = ape; pmse2 = ape;

for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, pred2{i}, predInt{i}, predInt2{i}, prob{i}, prob2{i}, R{i},...
        R2{i}, RInt{i}, RInt2{i}] = getNBAPETwoStep(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian 1 step
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end-1) - pred{i}).^2);
    
    % APE and predictive MSE for Bayesian 2 step
    ape2(i) = -sum(log(prob2{i}));
    pmse2(i) = mean((Iday(3:end) - pred2{i}).^2);
    
    disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
end

% Combine metrics, get best models
metrics = [ape' ape2' pmse' pmse2'];
[mins, mods] = min(metrics);
% Best windows sizes 
kbest = ks(mods);
disp(['Bayes k = [' num2str(kbest) ']' ]);


%% Incidence predictions and R estimates

% All metrics vs ks
figure;
subplot(2, 1, 1);
plot(ks, ape, 'color', grey2, 'linewidth', 2);
hold on; h = gca;
plot(ks, ape2, 'color', 'c', 'linewidth', 2);
plot([kbest(1) kbest(1)], h.YLim, '--', 'color', grey2, 'linewidth', 2);
plot([kbest(2) kbest(2)], h.YLim, '--', 'color', 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('APE');
xlim([ks(1) ks(end)]);
subplot(2, 1, 2);
plot(ks, pmse, 'color', grey2, 'linewidth', 2);
hold on; h = gca;
plot(ks, pmse2, 'color', 'c', 'linewidth', 2);
plot([kbest(3) kbest(3)], h.YLim, '--', 'color', grey2, 'linewidth', 2);
plot([kbest(4) kbest(4)], h.YLim, '--', 'color', 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('PMSE');
xlim([ks(1) ks(end)]);
xlabel('$k$ (days)');

% Optimal R windows and I predictions
Ropt = R{mods(1)}; Ropt2 = R{mods(2)};
dRopt = RInt{mods(1)}; dRopt2 = RInt{mods(2)};
Iopt = pred{mods(1)}; Iopt2 = pred{mods(2)};
dIopt = predInt{mods(1)}; dIopt2 = predInt{mods(2)};

% Plot from second time onwards
tplt = tday(2:end-1); tplt2 = tday(3:end);
Rplt = Rtrue(2:end-1); Rplt2 = Rtrue(3:end);
Iplt = Iday(2:end-1); Iplt2 = Iday(3:end);


% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
e1 = Ropt - dRopt(1, :); e1 = e1';
e2 = dRopt(2, :) - Ropt; e2 = e2';
plotCI(tplt', Ropt', e1, e2, 'c');
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ylim([0 4]);
subplot(2, 1, 2);
e1 = Iopt - dIopt(1, :); e1 = e1';
e2 = dIopt(2, :) - Iopt; e2 = e2';
plotCI(tplt', Iopt', e1, e2, 'c');
hold on;
%plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', 20);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('time (days)');
ylabel(['$k^* = $ ' num2str(kbest(1))]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['pred1step_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
e1 = Ropt2 - dRopt2(1, :); e1 = e1';
e2 = dRopt2(2, :) - Ropt2; e2 = e2';
plotCI(tplt2', Ropt2', e1, e2, 'c');
hold on;
plot(tplt2, Rplt2, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(2))]);
xlim([tplt2(1) tplt2(end)]);
ylim([0 4]);
subplot(2, 1, 2);
e1 = Iopt2 - dIopt2(1, :); e1 = e1';
e2 = dIopt2(2, :) - Iopt2; e2 = e2';
plotCI(tplt2', Iopt2', e1, e2, 'c');
hold on;
%plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', 20);
scatter(tplt2, Iplt2, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt2(1) tplt2(end)]);
xlabel('time (days)');
ylabel(['$k^* = $ ' num2str(kbest(2))]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['pred2step_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% Compare  mean predictions
figure;
e1 = Iopt(2:end) - dIopt(1, 2:end); e1 = e1';
e2 = dIopt(2, 2:end) - Iopt(2:end); e2 = e2';
plotCI(tplt(2:end)', Iopt(2:end)', e1, e2, 'b');
hold on;
e1 = Iopt2(1:end-1) - dIopt2(1, 1:end-1); e1 = e1';
e2 = dIopt2(2, 1:end-1) - Iopt2(1:end-1); e2 = e2';
plotCI(tplt(2:end)', Iopt2(1:end-1)', e1, e2, 'r');
grid off; box off; hold off;
xlim([tplt2(1) tplt2(end-1)]);
xlabel('time (days)');

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save([scenChoice '_' num2str(nks) '.mat']);
    cd(thisDir);
end


