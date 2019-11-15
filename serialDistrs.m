function pdistr = serialDistrs(tmax, distType)

% Assumptions and notes
% Chooses between discrete distributions on days
% Insert max days and calculate serial probabilities
% p must be a parameter, tday an array of integers


switch(distType)
    case 1
        % Geometric distribution over tmax for a given p
        pdistr = @(p) geomDistr(p, 1:tmax);
    case 2
        % Gamma distribution with integer shape (Erlang)
        pdistr = @(p) gammaDistr(p, 1:tmax);
    otherwise
        disp('No valid distribution specified');
        return;
end


% Gamma distribution, p is 1/mean, shape param preset
function pr = gammaDistr(p, k)

% Shape parameter
shapePm = 20;
% Scale parameter based on mean = shapePm*scalePm
scalePm = 1/(p*shapePm); ratePm = 1/scalePm;

% Gamma (Erlang) probabilities
pr = -log(factorial(shapePm-1)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(k) - ratePm*k;
pr = exp(pr);
%pr = gampdf(k, shapePm, scalePm);


% Geometric distribution, p is prob success
function pr = geomDistr(p, k)

% Set p to same dimension as k
p1 = p;
p = p*ones(size(k));

% Geometric distribution starting at 1 (vs 0)
pr = realpow((1-p), (k-1));
pr = p1*pr;