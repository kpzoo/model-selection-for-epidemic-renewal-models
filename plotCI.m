% Plot confidence intervals (asymmetric)
function [l, p] = plotCI(t, y, e1, e2, colstr)

% Assumptions and notes
% - uses boundedline package (assumed on path)
% - all inputs are column vectors, colstr = 'b' for example

% Main CI and median (or mean) plot
[l,p] = boundedline(t, y, [e1 e2], ['-' colstr]);

% Median line width and contrast of CI
l.LineWidth = 2; p.FaceAlpha = 0.5;
outlinebounds(l,p);