% Given look-back k, run one step ahead Bayesian predictions
function [pred, prob, R, wins] = getNegBinSwQuick(k, tPts, I, Lam)

% Assumptions and notes
% - vectorised for speed, removed interval computations
% - posterior for R based on k days look-back
% - if under k values available just use what is available

% Gamma prior parameters (a = shape, b = scale)
a = 1; b = 5; % from Cori 2013
%a = 10^(-5); b = 1000; % from Cauchemez 2006b
% Range of index time points
ir = 1:tPts-1; lenir = length(ir);

% Grouped incidence and infectiousness
B = zeros(1, lenir); A = B;

% At each time (tPts) compute historical R, predict at tPts+1
for i = ir
    % Look-back window of k (or less)
    idback = i:-1:max(i-k, 1); 
    % Relevant incidence sum (B) and total infectiousness sum (A)
    B(i) = sum(I(idback)); A(i) = sum(Lam(idback));
end

% Take length of last window (likely to be complete)
wins = length(idback);

% Parameters of posterior gamma on R
num = a + B;
den = 1/b + A;

% Posterior mean of R 
R = num./den;

% Predictive point for R (mean)
pred = Lam(ir).*R; % Lam(i) = Lam_{t+1}

% Prob of next incidence value - Pilatowska
p = Lam(ir)./den; p = p./(p + 1);
prob = nbinpdf(I(ir+1), num, 1-p);

% Check for numerical errors
if any(prob) == 0 || max(R) > 20 || any(isinf(prob))
    assignin('base', 'A', A);
    assignin('base', 'B', B);
    assignin('base', 'R', R);
    assignin('base', 'prob', prob);
    error('Issue with A or B or prob');
end






