% Combine the logSum values from 4 trajectories
clearvars; clc;
close all; 

% Assumption and notes
% - data must be in saveFol and cover all scenNam
% - data must have name matching scenNam

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Possible scenarios
scenNam = {'explosion', 'control', 'recovery', 'cascade'};
nScen = length(scenNam);

% Folder with datasets
saveFol = 'switch data'; thisDir = cd;
% Save figure
saveTrue = 1;

% Create a figure and plot logSum for each dataset
figure;
for i = 1:nScen
    % Data panel
    subplot(ceil(nScen/2), 2, i);
    % Load data of interest
    cd(saveFol);
    file = dir([scenNam{i} '*']);
    % Ensure unqique
    if length(file) ~= 1
        error('Too many files with same identifier');
    else
        load(file.name, 'treal', 'logSum', 'M', 'ts');
    end
    cd(thisDir);
    
    % Plot logSum over time and k
    plot(treal, logSum(:, 2:2:end-1), 'Color', grey1, 'LineWidth', 2);
    hold on;
    % Largest and smallest k case
    plot(treal, logSum(:, 1), 'Color', 'r', 'LineWidth', 2);
    plot(treal, logSum(:, end), 'Color', 'c', 'LineWidth', 2);
    % Switch times
    nSw = length(ts);
    h = gca; yl = h.YLim;
    for j = 1:nSw
        plot([ts(j) ts(j)], yl, '--', 'Color', grey2, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(scenNam{i});
    if i > nScen - 2
        xlabel('$t$ (days)');
    end
    xlim([treal(1) treal(end)]);
end

% Save combined figure
if saveTrue
    cd(saveFol);
    saveas(gcf, ['logSumsComb_' num2str(M)], 'fig');
    cd(thisDir);
end