% APE for prediction and model selection in renewal models
clearvars; clc;
close all; tic;

% Assumptions and notes
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
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};
scenNo = 5; scenChoice = scenNam{scenNo};

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
pred = cell(1, nks); predInt = pred;
% Posterior estimates of R over ks
R = pred; RInt = pred;
% APE metric and PMSE
prob = pred; ape = zeros(1, nks); pmse = ape;
% Parameters of NB predictive distribution
pm1 = pred; pm2 = pred;

for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, pm1{i}, pm2{i}] = getNegBinAPEData(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
    
    disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
end

% Best models according to metrics
[apeMin, apeMod] = min(ape);
[pmseMin, pmseMod] = min(pmse);

% Best ks and nGrps
modID = [apeMod pmseMod];
kbest = ks(modID);
disp(['Bayes k: [ape pmse] = [' num2str(kbest) ']' ]);

% Accuracy of predictions over ks
percMiss = zeros(1, nks);
Itrue = Iday(2:end);
for i = 1:nks
    predI = predInt{i};
    % Ids at which incidence outside credible interval
    idout = union(find(Itrue < predI(1, :)), find(Itrue > predI(2, :)));
    % Percentage missed
    percMiss(i) = 100*length(idout)/length(Itrue);
end

%% Incidence predictions, R estimates, main figures

% For plotting set k to window length = k + 1
ks = ks + 1; kbest = kbest + 1;

% Overfit, underfit and optimal R windows
Rov = R{1}; Rund = R{end}; Ropt = R{modID(1)};
% Plot from second time onwards
tplt = tday(2:end); Rplt = Rtrue(2:end);
% Best NB predictive parameters
pmbest1 = pm1{apeMod}; pmbest2 = pm2{apeMod};

if testFig

    % Model selection with APE and PMSE
    figure;
    yyaxis left
    plot(ks, ape,'linewidth', 2);
    hold on;
    plot(kbest(1), apeMin, 'o', 'MarkerSize', 10);
    hold off;
    ylabel('APE');
    yyaxis right
    plot(ks, pmse, 'linewidth', 2);
    hold on;
    plot(kbest(2), pmseMin, 'o', 'MarkerSize', 10);
    hold off;
    ylabel('PMSE');
    grid off; box off;
    xlabel('$k$ (days)');
    
    % PMSE against predictive error
    figure;
    yyaxis left
    plot(ks, pmse, '.-', 'linewidth', 2, 'markersize', 20);
    h = gca; h.YColor = h.XColor;
    yyaxis right
    plot(ks, percMiss, '.-', 'linewidth', 2, 'markersize', 20);
    h = gca; h.YColor = h.XColor;
    grid off; box off;
    legend('PMSE', 'predictive error', 'location', 'best');
    xlabel('$k$ (days)');
    
    % One step ahead predictions and metric (APE vs PMSE)
    figure;
    subplot(2, 2, 1);
    e1 = Ropt - RInt{apeMod}(1, :); e1 = e1';
    e2 = RInt{apeMod}(2, :) - Ropt; e2 = e2';
    plotCI(tplt', Ropt', e1, e2, 'c');
    hold on;
    plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
    grid off; box off; hold off;
    ylabel(['$k^* = $ ' num2str(kbest(1))]);
    xlim([tplt(1) tplt(end)]);
    ylim([0 4]);
    subplot(2, 2, 3);
    e1 = pred{apeMod} - predInt{apeMod}(1, :); e1 = e1';
    e2 = predInt{apeMod}(2, :) - pred{apeMod}; e2 = e2';
    plotCI(tplt', pred{apeMod}', e1, e2, 'c');
    hold on;
    %plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', 20);
    scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    grid off; box off; hold off;
    xlim([tplt(1) tplt(end)]);
    xlabel('time (days)');
    ylabel(['$k^* = $ ' num2str(kbest(1))]);
    
    subplot(2, 2, 2);
    e1 = Ropt - RInt{pmseMod}(1, :); e1 = e1';
    e2 = RInt{pmseMod}(2, :) - Ropt; e2 = e2';
    plotCI(tplt', Ropt', e1, e2, 'c');
    hold on;
    plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
    grid off; box off; hold off;
    ylabel(['$k^* = $ ' num2str(kbest(2))]);
    xlim([tplt(1) tplt(end)]);
    ylim([0 4]);
    subplot(2, 2, 4);
    e1 = pred{pmseMod} - predInt{pmseMod}(1, :); e1 = e1';
    e2 = predInt{pmseMod}(2, :) - pred{pmseMod}; e2 = e2';
    plotCI(tplt', pred{pmseMod}', e1, e2, 'c');
    hold on;
    %plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', 20);
    scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    grid off; box off; hold off;
    xlim([tplt(1) tplt(end)]);
    xlabel('time (days)');
    ylabel(['$k^* = $ ' num2str(kbest(2))]);
    if saveTrue
        cd(saveFol);
        saveas(gcf, ['PMSEAPE_' scenChoice '_' num2str(nks)], 'fig');
        cd(thisDir);
    end
    
    % Indices of interest
    lenx = 20; x = 1:2*max(Iday);
    idx = unique(round(linspace(1, length(pmbest1), lenx)));
    lenx = length(idx); nb = cell(1, lenx);
    
    % All NB predictive distributions
    figure;
    hold on;
    for i = 1:lenx
        % Centre relative to mean
        [xmean, xvar] = nbinstat(pmbest1(idx(i)), pmbest2(idx(i)));
        % Domain to consider
        %x = 1:xmean+3*sqrt(xvar);
        % NB values
        nb{i} = nbinpdf(x, pmbest1(idx(i)), pmbest2(idx(i)));
        plot(x-xmean, nb{i}, 'linewidth', 2);
    end
    grid off; box off; hold off;
    xlabel('$I_{t+1}$');
    ylabel('$P(I_{t+1})$');
    if saveTrue
        cd(saveFol);
        saveas(gcf, ['distr_' scenChoice '_' num2str(nks)], 'fig');
        cd(thisDir);
    end
    
    % Parameters across time
    figure;
    subplot(2, 1, 1);
    plot(tplt, pmbest1, 'linewidth', 2);
    ylabel('$r$');
    grid off; box off; 
    subplot(2, 1, 2);
    plot(tplt, pmbest2, 'linewidth', 2);
    ylabel('$p$');
    grid off; box off;
    xlabel('time (days)');
end

% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
e1 = Ropt - RInt{apeMod}(1, :); e1 = e1';
e2 = RInt{apeMod}(2, :) - Ropt; e2 = e2';
plotCI(tplt', Ropt', e1, e2, 'c');
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ylim([0 4]);
subplot(2, 1, 2);
e1 = pred{apeMod} - predInt{apeMod}(1, :); e1 = e1';
e2 = predInt{apeMod}(2, :) - pred{apeMod}; e2 = e2';
plotCI(tplt', pred{apeMod}', e1, e2, 'c');
hold on;
%plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', 20);
scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('time (days)');
ylabel(['$k^* = $ ' num2str(kbest(1))]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['predRhat_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

figure;
subplot(3, 2, 1);
e1 = Ropt - RInt{apeMod}(1, :); e1 = e1';
e2 = RInt{apeMod}(2, :) - Ropt; e2 = e2';
plotCI(tplt', Ropt', e1, e2, 'c');
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ylim([0 4]);
subplot(3, 2, 3);
e1 = Rund - RInt{end}(1, :); e1 = e1';
e2 = RInt{end}(2, :) - Rund; e2 = e2';
plotCI(tplt', Rund', e1, e2, 'c');
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(end))]);
xlim([tplt(1) tplt(end)]);
ylim([0 4]);
subplot(3, 2, 5);
e1 = Rov - RInt{1}(1, :); e1 = e1';
e2 = RInt{1}(2, :) - Rov; e2 = e2';
plotCI(tplt', Rov', e1, e2, 'c');
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(1))]);
xlim([tplt(1) tplt(end)]);
ylim([0 4]);
xlabel('time (days)');

mksz = 10;
ax(1) = subplot(3, 2, 2);
e1 = pred{apeMod} - predInt{apeMod}(1, :); e1 = e1';
e2 = predInt{apeMod}(2, :) - pred{apeMod}; e2 = e2';
plotCI(tplt', pred{apeMod}', e1, e2, 'c');
hold on;
plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', mksz);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(3, 2, 4);
e1 = pred{end} - predInt{end}(1, :); e1 = e1';
e2 = predInt{end}(2, :) - pred{end}; e2 = e2';
plotCI(tplt', pred{end}', e1, e2, 'c');
hold on;
plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', mksz);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(end))]);
xlim([tplt(1) tplt(end)]);
ax(3) = subplot(3, 2, 6);
e1 = pred{1} - predInt{1}(1, :); e1 = e1';
e2 = predInt{1}(2, :) - pred{1}; e2 = e2';
plotCI(tplt', pred{1}', e1, e2, 'c');
hold on;
plot(tplt, Iday(2:end), '.', 'color', grey1, 'linewidth', 2, 'markersize', mksz);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(1))]);
xlim([tplt(1) tplt(end)]);
xlabel('time (days)');
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['compRI_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% APE against predictive error
figure;
yyaxis left
plot(ks, ape, '.', 'linewidth', 2, 'markersize', 20);
hold on; h = gca;
plot([kbest(1) kbest(1)], h.YLim, 'k--', 'LineWidth', 2);
hold off; 
h = gca; h.YColor = h.XColor;
yyaxis right
plot(ks, percMiss, '.', 'linewidth', 2, 'markersize', 20);
h = gca; h.YColor = h.XColor;
grid off; box off;
%legend('APE', 'predictive error', 'location', 'best');
xlabel('$k$ (days)');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['metric_' scenChoice '_' num2str(nks)], 'fig');
    cd(thisDir);
end


% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save([scenChoice '_' num2str(nks) '.mat']);
    cd(thisDir);
end


