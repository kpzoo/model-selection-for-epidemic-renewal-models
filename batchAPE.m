% Run many epidemics and APE forecasts
clearvars; clc;
close all; tic;

% Assumptions and notes
% - uses NB posterior predictions and Bayesian APE
% - APE compares sequential predictions with true values
% - simulates a single epidemic, predicts I(t), estimates R(t)

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving
saveFol = 'batch data'; thisDir = cd;

% Define a scenario
scenNam = {'constant', 'cyclic', 'logistic', 'piecewise', 'boom-bust', 'bottle', '2-step', 'filtered'};
scenNo = 2; 
scenChoice = scenNam{scenNo};

% Num of independent epidemics
M = 100;
% Starting times (days) 
tday0 = 1:201; nday0 = length(tday0);

% Window lengths to search across
ks = 1:ceil(nday0/2); nks = length(ks);
disp(['k varies from ' num2str(ks(1)) ' to ' num2str(ks(end))]);

% Outputs
Rmse = cell(1, M); pMiss = Rmse; 
modID = zeros(1, M); ape = Rmse; Imse = Rmse; 

% Simulate epidemics and select kbest
for i = 1:M
    % Main APE forecasting
    [pMiss{i}, modID(i), ape{i}, Rmse{i}, Imse{i}] = apeForecastFn(scenNo,...
        ks, nks, tday0, nday0);
    disp(['Completed: ' num2str(i) ' of ' num2str(M)]);
end

% Best windows
kbest = ks(modID); kavg = mean(kbest);
disp(['Mean k^* = ' num2str(kavg)]);

% Convert to matrix
Rmse = cell2mat(Rmse'); Imse = cell2mat(Imse');
pMiss = cell2mat(pMiss'); ape = cell2mat(ape');

% Average ape, Imse etc over runs
apeMed = quantile(ape, 0.5); apeInt = quantile(ape, [0.0275 0.975]);
percMed = quantile(pMiss, 0.5); percInt = quantile(pMiss, [0.0275 0.975]);
dIMed = quantile(Imse, 0.5); dIInt = quantile(Imse, [0.0275 0.975]);
dRMed = quantile(Rmse, 0.5); dRInt = quantile(Rmse, [0.0275 0.975]);

% Best k for other metrics (medians)
ids = zeros(1, 4);
[~, ids(1)] = min(apeMed); [~, ids(2)] = min(percMed);
[~, ids(3)] = min(dIMed); [~, ids(4)] = min(dRMed);
kmets = ks(ids);

% All APE windows histogram counts
figure;
% Set so recover ks directly in terms of edges
[kcounts, kvals] = histcounts(kbest, [ks ks(end)+1]);
kvals = kvals(1:end-1);
h = stem(kvals, kcounts);
h.Marker = '.'; h.MarkerSize = 30; h.Color = 'c';
h.LineWidth = 2; h.MarkerFaceColor = 'c'; h.MarkerEdgeColor = 'c';
xlim([ks(1) ks(end)]);
grid off; box off;
xlabel('$k$'); ylabel('freq');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['kbest_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% Plot performance measures
figure;
subplot(2, 2, 1);
plot(ks, apeMed, 'color', grey1, 'linewidth', 2);
hold on; xlim([ks(1) ks(end)]);
plot(ks, apeInt, 'color', grey2, 'linewidth', 2);
h = gca; plot([kmets(1) kmets(1)], h.YLim, 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('APE');

subplot(2, 2, 2);
plot(ks, dIMed, 'color', grey1, 'linewidth', 2);
hold on; xlim([ks(1) ks(end)]);
plot(ks, dIInt, 'color', grey2, 'linewidth', 2);
h = gca; plot([kmets(1) kmets(1)], h.YLim, 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('mse(I)');

subplot(2, 2, 3);
plot(ks, dRMed, 'color', grey1, 'linewidth', 2);
hold on; xlim([ks(1) ks(end)]);
plot(ks, dRInt, 'color', grey2, 'linewidth', 2);
h = gca; plot([kmets(1) kmets(1)], h.YLim, 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('mse(R)'); xlabel('$k$ (days)');

subplot(2, 2, 4);
plot(ks, percMed, 'color', grey1, 'linewidth', 2);
hold on; xlim([ks(1) ks(end)]);
plot(ks, percInt, 'color', grey2, 'linewidth', 2);
h = gca; plot([kmets(1) kmets(1)], h.YLim, 'c', 'linewidth', 2);
hold off; grid off; box off;
ylabel('Miss \%'); xlabel('$k$ (days)');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['allMets_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save([scenChoice '_' num2str(M) '_' num2str(nks) '.mat']);
    cd(thisDir);
end



