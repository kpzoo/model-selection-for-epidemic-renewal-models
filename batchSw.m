% Run many epidemics and APE forecasts
clearvars; clc;
close all; tic;

% Assumptions and notes
% - gets successive k* choices across time in batch
% - keeps mean estimates and predictions

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

%% Setup simulations

% Choose a scenario
scenNo = 4;

% Save data and test figs
saveTrue = 0;
% Folder for saving
saveFol = 'switch data'; thisDir = cd;

% Define all scenarios
scenNam = {'explosion', 'control', 'recovery', 'cascade'};
scenChoice = scenNam{scenNo};
switch(scenNo)
    % Rs are distinct reprod nums, ts are switch points
    case 1
        % Exploding epidemic
        Rs = [0.8 2];
        ts = 50;
    case 2
        % Rapidly controlled epidemic
        Rs = [2 0.8];
        ts = 100;
    case 3
        % Rapid control that recovers
        Rs = [2 0.8 1.2];
        ts = [80 110];
    case 4
        % Two stage control
        Rs = [5 1.5 0.8];
        ts = [40 120];
end

% Num of independent epidemics
M = 1000;
% Starting times (days)
tday0 = 1:151; nday0 = length(tday0);

% Window lengths to search across
ks = 1:2:99;
nks = length(ks);
disp(['k varies from ' num2str(ks(1)) ' to ' num2str(ks(end))]);

% Storage variables
s = cell(1, M); z = zeros(1, M);
% Successive ks ids for best models
apeIDs = s; pmseIDs = s;
% True trajectory for R and I
Rtrue = s; tday = s; Itrue = s;
% Metric values across ks
ape = s; pmse = s;
% Best metric id for each run
apei = z; pmsei = z;
% Mean R estimates at best metric
Rhat = s; Rhatpmse = s;
% Mean I predictions at best metric
Ihat = s; Ihatpmse = s;
% Start and end times
tstart = z; tend = z;
% Prediction quality (log and mse)
logPred = s; msePred = s;

%% Compute APE and other metrics

% Simulate epidemics and select kbest
for i = 1:M
    % Main APE forecasting
    [apeIDs{i}, pmseIDs{i}, tday{i}, Rhat{i}, Rhatpmse{i}, Ihat{i}, Ihatpmse{i}, ape{i},...
        pmse{i}, apei(i), pmsei(i), Rtrue{i}, Itrue{i}, logPred{i}, msePred{i}]...
        = apeForecastSwitchFn(ts, Rs, tday0, ks, nks, scenNo);
    disp(['Completed: ' num2str(i) ' of ' num2str(M)]);
    
    % Get common times of each run
    tstart(i) = tday{i}(1); tend(i) = tday{i}(end);
end

% Adjust for true window length
ks = ks + 1;
% Extract best metrics and ids
ape = cell2mat(ape'); pmse = cell2mat(pmse');
kape = ks(apei); kpmse = ks(pmsei);

% Range of time common to all epidemics
t = max(tstart):min(tend); lent = length(t);
% True R over this time
R = Rtrue{1}(ismember(tday{1}, t));
Rmin = min(R); Rmax = max(R); clc;
disp(['R between ' [num2str(Rmin) ' ' num2str(Rmax)]]);
% Exclude last t value in prediction
treal = t(1):t(end-1);

% Store truncated variables in this matrix
zM = zeros(lent-1, M);
% Actual truncated declarations
apeIDsM = zM; pmseIDsM = zM; ItrueM = zM;
RhatM = zM; RhatM2 = zM; IhatM = zM; IhatM2 = zM;

% Get indices in tday corresponding to these and truncate
for i = 1:M
    % Valid range
    v = ismember(tday{i}, treal); v = find(v);
    % Truncate variables
    apeIDsM(:, i) = apeIDs{i}(v); pmseIDsM(:, i) = pmseIDs{i}(v);
    RhatM(:, i) = Rhat{i}(v); RhatM2(:, i) = Rhatpmse{i}(v);
    IhatM(:, i) = Ihat{i}(v); IhatM2(:, i) = Ihatpmse{i}(v);
    % For true incidence only look at values to be predicted
    ItrueM(:, i) = Itrue{i}(v);
    % Fix time range of predictive scores
    logPred{i} = logPred{i}(v, :);
    msePred{i} = msePred{i}(v, :);
end

% Window sizes corresponding to APE and PMSE
kapeIDsM = ks(apeIDsM); kpmseIDsM = ks(pmseIDsM);

% Get quantiles of estimates and predictions
kapeM = quantile(kapeIDsM', [0.025, 0.5, 1-0.025]);
kpmseM = quantile(kpmseIDsM', [0.025, 0.5, 1-0.025]);
RvM = quantile(RhatM', [0.025, 0.5, 1-0.025]);
RvM2 = quantile(RhatM2', [0.025, 0.5, 1-0.025]);
% Use complete range for incidence
IvM = quantile(IhatM', [0, 0.5, 1]);
IvM2 = quantile(IhatM2', [0, 0.5, 1]);

% Overall choice of k*
kst_ape = kapeM(2, end); kst_pmse = kpmseM(2, end);
disp(['k*(end) for [ape pmse] = ' [num2str(kst_ape) ' ' num2str(kst_pmse)]]);

% Get size of matrices in cell.
matSize = size(logPred{1},1);
% Create 3D matrix for scores
logA = reshape(cell2mat(logPred), matSize,[], M);
mseA = reshape(cell2mat(msePred), matSize,[], M);
% Sum 3D matrix along 3rd dimension
logSum = sum(logA, 3); mseSum = sum(mseA, 3);

%% Save data and figures

% Optimal k* by both methods
figure;
h = histogram(kape, ks);
hold on;
histogram(kpmse, ks);
hold off; grid off; box off;
xlabel('$k$ (days)');
ylabel('freq($k^*$)');
legend('APE', 'PMSE', 'location', 'best');

% Predictions in total
figure;
plot(treal, ItrueM, '.', 'Color', grey1, 'MarkerSize', 20);
hold on;
[~, em] = stairs(treal, IvM(2, :));
e1 = IvM(2, :) - IvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = IvM(3, :) - IvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
grid off; box off;
ylabel('$I_t$');
xlabel('$t$ (days)');
xlim([treal(1) treal(end)]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['apePred_' num2str(M) '_' num2str(scenNo)], 'fig');
    cd(thisDir);
end

% Predictive scroes as a function of ks and time
figure;
plot(treal, logSum, 'Color', grey1, 'LineWidth', 2);
hold on;
plot(treal, logSum(:, 1), 'Color', 'r', 'LineWidth', 2);
plot(treal, logSum(:, end), 'Color', 'm', 'LineWidth', 2);
plot(treal, mseSum, 'Color', grey1, 'LineWidth', 2);
plot(treal, mseSum(:, 1), 'Color', 'r', 'LineWidth', 2);
plot(treal, mseSum(:, end), 'Color', 'm', 'LineWidth', 2);
hold off; grid off; box off;
ylabel('predictive score');
xlabel('$t$ (days)');
xlim([treal(1) treal(end)]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['scorePred_' num2str(M) '_' num2str(scenNo)], 'fig');
    cd(thisDir);
end

% Confidence across time for pmse
figure;
% R estimates with optimal grouping
h(1) = subplot(2, 1, 1);
[~, em] = stairs(treal, RvM(2, :));
e1 = RvM(2, :) - RvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = RvM(3, :) - RvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
hold on;
stairs(treal, R(1:end-1), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
ylabel('$R_{\tau}$');
xlim([treal(1) treal(end)]);
% Successive k* choices until optimal at end
h(2) = subplot(2, 1, 2);
[~, em] = stairs(treal, kapeM(2, :));
e1 = kapeM(2, :) - kapeM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = kapeM(3, :) - kapeM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
grid off; box off;
ylabel('$k^*$');
xlabel('$t$ (days)');
xlim([treal(1) treal(end)]);
ylim([ks(1) ks(end)+1]);
linkaxes(h, 'x');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['ape_' num2str(M) '_' num2str(scenNo)], 'fig');
    cd(thisDir);
end


figure;
% Incidence predictions and true from ape
h(1) = subplot(3, 1, 1);
plot(treal, ItrueM, '.', 'Color', grey1, 'MarkerSize', 20);
hold on;
[~, em] = stairs(treal, IvM(2, :));
e1 = IvM(2, :) - IvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = IvM(3, :) - IvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
grid off; box off;
ylabel('$I_t$');
xlim([treal(1) treal(end)]);
h(2) = subplot(3, 1, 2);
% R estimates with APE grouping
[~, em] = stairs(treal, RvM(2, :));
e1 = RvM(2, :) - RvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = RvM(3, :) - RvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
hold on;
% True R
stairs(treal, R(1:end-1), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
ylabel('$R_{\tau}$');
xlim([treal(1) treal(end)]);
% Successive k* choices from APE
h(3) = subplot(3, 1, 3);
[~, em] = stairs(treal, kapeM(2, :));
e1 = kapeM(2, :) - kapeM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = kapeM(3, :) - kapeM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
hold off; grid off; box off;
ylabel('$k^*$');
xlabel('$t$ (days)');
xlim([treal(1) treal(end)]);
ylim([ks(1) ks(end)+1]);
linkaxes(h, 'x');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['comb_' num2str(M) '_' num2str(scenNo)], 'fig');
    cd(thisDir);
end

% All plots combined
figure;
% Incidence predictions and true from ape
h(1) = subplot(3, 1, 1);
plot(treal, ItrueM, '.', 'Color', grey1, 'MarkerSize', 20);
hold on;
[~, em] = stairs(treal, IvM(2, :));
e1 = IvM(2, :) - IvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = IvM(3, :) - IvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
grid off; box off;
ylabel('$I_t$');
xlim([treal(1) treal(end)]);
% R estimates with PMSE grouping
h(2) = subplot(3, 1, 2);
[~, em] = stairs(treal, RvM2(2, :));
e1 = RvM2(2, :) - RvM2(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = RvM2(3, :) - RvM2(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
hold on;
% R estimates with APE grouping
[~, em] = stairs(treal, RvM(2, :));
e1 = RvM(2, :) - RvM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = RvM(3, :) - RvM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'r');
% True R
stairs(treal, R(1:end-1), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
ylabel('$R_{\tau}$');
xlim([treal(1) treal(end)]);
% Successive k* choices from PMSE
h(3) = subplot(3, 1, 3);
[~, em] = stairs(treal, kpmseM(2, :));
e1 = kpmseM(2, :) - kpmseM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = kpmseM(3, :) - kpmseM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'c');
hold on;
% Successive k* choices from APE
[~, em] = stairs(treal, kapeM(2, :));
e1 = kapeM(2, :) - kapeM(1, :); e1 = e1';
[~, e1] = stairs(treal, e1);
e2 = kapeM(3, :) - kapeM(2, :); e2 = e2';
[tr, e2] = stairs(treal, e2);
plotCI(tr, em, e1, e2, 'r');
hold off; grid off; box off;
ylabel('$k^*$');
xlabel('$t$ (days)');
xlim([treal(1) treal(end)]);
ylim([ks(1) ks(end)+1]);
linkaxes(h, 'x');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['comball_' num2str(M) '_' num2str(scenNo)], 'fig');
    cd(thisDir);
end


% Timing
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Save data
if saveTrue
    cd(saveFol);
    % List of variables to save
    listVar = {'treal', 'kapeM', 'kpmseM', 'ks', 'R', 'Rs', 'ts', 'scenNo', 'scenChoice',...
        'RvM', 'RvM2', 'IvM', 'IvM2', 'ItrueM', 'mseSum', 'logSum', 'kpmse', 'kape', 'M', 'nks'};
    save([scenChoice '_' num2str(M) '_' num2str(scenNo) '.mat'], listVar{:});
    cd(thisDir);
end

