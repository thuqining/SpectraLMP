function pol = finish_policy(pol, xtrain, par, beta, fit)
% FINISH_POLICY  Algorithm 1 (v23) from a synthesized FIR bank:
%   - gain beta on the taps;
%   - no filter output for units whose band lies below the resolvable
%     frequency 2*pi/(L*dt) rad/h of an L-tap filter (their upper band edge
%     is the crossover with the next faster unit);
%   - stored-energy values for step 3, fitted on the same training window
%     (fit_values), unless a fit is supplied.
pol.theta = beta * pol.theta;
cx = sqrt(max(diff(-pol.kappa), 0)./max(diff(pol.ceff), eps));    % rad/h
pol.unresolved = find([inf, cx] < 2*pi/(pol.L*par.dt));
pol.theta(:, pol.unresolved) = 0;
if nargin < 5 || isempty(fit)
    fit = fit_values(xtrain, par);
end
pol.valueFit = fit;
end
