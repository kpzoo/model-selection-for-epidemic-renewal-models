% Given look-back k, run one step ahead Bayesian predictions
function [pred, pred2, predInt, predInt2, prob, prob2, R, R2, RInt, RInt2] = ...
    getNBAPETwoStep(k, tPts, I, Lam)

% Assumptions and notes
% - updated for iterated 2 step ahead forecasts
% - posterior for R based on k days look-back
% - if under k values available just use what is available
% - point and interval prediction from negative binomial 
% - conjugate gamma priors used for renewal model

% Gamma prior parameters (a = shape, b = scale)
a = 1; b = 5; % from Cori 2013

% Range of index time points (shortened to get t+2)
ir = 1:tPts-2; imax = length(ir);

% Grouped incidence and infectiousness
B = zeros(1, imax); A = B;
% Confidence intervals
predInt = zeros(2, imax); RInt = predInt;

% At each time (tPts) compute historical R, predict at tPts+1
for i = ir
    % Look-back window of k (or less) - note to i-k vs i-k+1
    idback = i:-1:max(i-k, 1);
    % Relevant incidence sum (B) and total infectiousness sum (A)
    B(i) = sum(I(idback)); A(i) = sum(Lam(idback));
end

% Parameters of posterior gamma on R
num = a + B;
den = 1/b + A;

% Posterior mean of R (1 step ahead)
R = num./den;
% Iterated R estimates
RInt(1, :) = gaminv(0.025, num, 1./den);
RInt(2, :) = gaminv(0.975, num, 1./den);

% Predictive 1 step ahead incidence
pred = Lam(ir).*R; % Lam(i) = Lam_{t+1}

% Negative binomial 1 step confidence
p = Lam(ir)./den; p = p./(p + 1);
predInt(1, :) = nbininv(0.025, num, 1-p);
predInt(2, :) = nbininv(0.975, num, 1-p);

% Prob of next incidence value - Pilatowska
prob = nbinpdf(I(ir+1), num, 1-p);

% Update the incidence part of R gamma
num2 = num + pred;
% Update infectiousness, Lam_{t+1} = f(I_t)
den2 = den + Lam(2:end-1);

% Posterior mean of R (2 step ahead)
R2 = num2./den2;
% Iterated R estimates
RInt2(1, :) = gaminv(0.025, num2, 1./den2);
RInt2(2, :) = gaminv(0.975, num2, 1./den2);

% Predictive 2 step ahead incidence
pred2 = Lam(ir+1).*R2; % Lam(i) = Lam_{t+1}

% Negative binomial 1 step confidence
p2 = Lam(ir+1)./den2; p2 = p2./(p2 + 1);
predInt2(1, :) = nbininv(0.025, num2, 1-p2);
predInt2(2, :) = nbininv(0.975, num2, 1-p2);

% Prob of 2 step incidence value 
prob2 = nbinpdf(I(ir+2), num, 1-p);





