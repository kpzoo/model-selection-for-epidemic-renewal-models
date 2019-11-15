% Combined metric plot

figure; 
%%
% Change subplot as construct

subplot(3, 2, 6); % chaange plot with each data


% Simply input mat file and copy this code

% APE against predictive error
yyaxis left
plot(ks, ape, '.', 'linewidth', 2, 'markersize', 10);
hold on; h = gca;
plot([kbest(1) kbest(1)], h.YLim, 'k--', 'LineWidth', 2);
hold off; ylabel(scenChoice);
h = gca; h.YColor = h.XColor; h.YTickLabel = '';
yyaxis right
plot(ks, percMiss, '.', 'linewidth', 2, 'markersize', 10);
h = gca; h.YColor = h.XColor; h.YTickLabel = '';
grid off; box off;
%legend('APE', 'predictive error', 'location', 'best');
xlabel('$k$ (days)');


%%
if saveTrue
    cd(saveFol);
    saveas(gcf, ['metricComb' '_' num2str(nks)], 'fig');
    cd(thisDir);
end