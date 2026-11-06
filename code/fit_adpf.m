function pa = fit_adpf(xtr, hrs, par, Ht, memH)
% FIT_ADPF  ADP benchmark with the same information as Algorithm 1: the
% marginal value of unit i is v_i(e, f_t) = max(theta_i' [1, f_t] + v1_i e, 0),
% v1_i <= 0, with f_t the standardized features of price_features (same
% memory memH and ridge penalty as design_price), fitted by least squares
% to the SoC duals of the training hindsight QP Ht against the SoC e_{i,t-1}.
n = numel(xtr); N = par.N;
F = price_features(xtr, hrs, memH);
idx = (memH*4+1):n;                            % full history available
mu = mean(F(idx,:)); sd = std(F(idx,:)) + eps;
Z = [ones(numel(idx),1), (F(idx,:) - mu)./sd];
k = size(Z, 2); lam = 0.1*numel(idx);
B = zeros(k, N); v1 = zeros(1, N);
for i = 1:N
    u = [par.e0(i), Ht.e(i,1:end-1)]'; u = u(idx);
    D = Ht.shadow_value(i,:)'; D = D(idx);
    A = [Z, u]; R = lam*diag([0, ones(1, k-1), 0]);
    th = (A'*A + R) \ (A'*D);
    if th(end) > 0                             % enforce v1 <= 0 (concave value)
        th = [(Z'*Z + R(1:k,1:k)) \ (Z'*D); 0];
    end
    B(:,i) = th(1:k); v1(i) = th(end);
end
pa = struct('B', B, 'v1', v1, 'mu', mu, 'sd', sd, 'memH', memH);
end
