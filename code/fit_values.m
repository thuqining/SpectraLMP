function fit = fit_values(xtrain, par, train)
% FIT_VALUES  Separable stored-energy values v_i(e) = max(alpha_i + beta_i e, 0),
% beta_i <= 0, fitted to the SoC duals of the training-window hindsight QP.
% This is the same fit as the ADP benchmark (method_adp); fit(:,i) holds
% [alpha_i; beta_i] in $/kWh and $/kWh^2. A precomputed
% method_hindsight(xtrain, par) may be passed as the third argument.
if nargin < 3, [~, train] = evalc('method_hindsight(xtrain, par)'); end
fit = zeros(2, par.N);
for i = 1:par.N
    u = [par.e0(i), train.e(i,1:end-1)]';
    v = train.shadow_value(i,:)';
    cu = u - mean(u);
    slope = min(sum(cu.*(v - mean(v)))/max(sum(cu.^2), eps), 0);
    fit(:,i) = [mean(v) - slope*mean(u); slope];
end
end
