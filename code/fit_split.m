function ps = fit_split(xtr, hrs, par, Ht, memH, tau, bus)
% FIT_SPLIT  Split predictor of the opportunity value (Theorem 1). The average
% (4) = causal half + anticausal half. The causal half is observed online as
% the exponential average Cm_i of the past LMPs at the unit's bus with time
% constant tau_i; the anticausal half is predicted from the past net load.
% Unit i: mu~_i = theta_i' [1, f_t, Cm_{i,t-1}], ridge-fitted (as design_price)
% to the SoC duals of the training hindsight QP Ht, with Cm from its LMPs.
n = numel(xtr); N = par.N; dt = par.dt;
F = price_features(xtr, hrs, memH);
idx = (memH*4+1):n;
mu = mean(F(idx,:)); sd = std(F(idx,:)) + eps;
Zf = [ones(numel(idx),1), (F(idx,:) - mu)./sd];
k = size(Zf, 2); lam = 0.1*numel(idx);
r = exp(-dt ./ tau(:)');
B = zeros(k, N); g = zeros(1, N); muC = zeros(1, N); sdC = ones(1, N); Cend = zeros(1, N);
for i = 1:N
    l = Ht.lmp(:, bus(i));
    Cm = filter(1 - r(i), [1, -r(i)], l, l(1)*r(i));     % average up to step t
    Cend(i) = Cm(end);
    Cl = [l(1); Cm(1:end-1)];                            % up to step t-1 (decision time)
    c = Cl(idx); muC(i) = mean(c); sdC(i) = std(c) + eps;
    A = [Zf, (c - muC(i))/sdC(i)];
    R = lam*diag([0, ones(1, k-1), 1]);
    th = (A'*A + R) \ (A'*Ht.shadow_value(i, idx)');
    B(:,i) = th(1:k); g(i) = th(end);
end
ps = struct('B', B, 'g', g, 'muC', muC, 'sdC', sdC, 'mu', mu, 'sd', sd, 'memH', memH, 'Cend', Cend, 'r', r);
end
