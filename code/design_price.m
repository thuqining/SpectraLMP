function pp = design_price(xtr, hrs, par, pol, opts)
% DESIGN_PRICE  Offline design of the price-form policy (Algorithm 1).
% Shape: each unit's cost average LP_i[lambda], LP_i(w) = kappa_i/(kappa_i
% + c_i w^2), of the training hindsight balancing cost (Theorem 1).
% Level (opts.lossCal): an affine map per unit, fitted to the SoC duals of
% the training QP, converts the cost average to the value of stored energy
% under losses and limits. Units whose time constant exceeds the predictor
% memory are valued at the window mean of their duals.
% Predictor: ridge regression on causal features (price_features), a
% finite-memory linear filter of the past net load plus daily harmonics.
%   opts.Ht        precomputed method_hindsight(xtr, par) (optional)
%   opts.noHarm    drop the daily harmonics (ablation)
%   opts.constant  predict the training mean only (ablation)
if nargin < 5, opts = struct(); end
dt = par.dt; N = par.N; n = numel(xtr);
m0 = 192; if isfield(opts, 'memH'), m0 = opts.memH; end   % default 8 days
memH = min(m0, floor(n*dt/2));               % memory <= half the training window
if isfield(opts, 'Ht'), Ht = opts.Ht; else, [~, Ht] = evalc('method_hindsight(xtr, par)'); end
f = (0:n-1)'/n; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
Y = zeros(n, N);
for i = 1:N
    pr = Ht.price(:);
    if isfield(Ht, 'priceUnit'), pr = Ht.priceUnit(:, i); end   % network: LMP at unit bus
    Y(:,i) = real(ifft(fft(pr) .* (pol.kappa(i)./(pol.kappa(i) + pol.ceff(i)*w.^2))));
end
useHarm = ~(isfield(opts, 'noHarm') && opts.noHarm);
cal = zeros(2, N); cal(2,:) = 1;
isDirect = isfield(opts, 'direct') && opts.direct;   % ablation: regress the hindsight SoC duals directly
if isDirect
    Y = Ht.shadow_value';                     % n x N, no averaging (4), no affine map, no constant rule
elseif isfield(opts, 'lossCal') && opts.lossCal
    tau = sqrt(pol.ceff./pol.kappa);          % h
    for i = 1:N
        D = Ht.shadow_value(i,:)';
        if tau(i) > memH                      % unresolvable: window mean
            cal(:,i) = [mean(D); 0];
        else
            cal(:,i) = [ones(n,1), Y(:,i)] \ D;
        end
        Y(:,i) = cal(1,i) + cal(2,i)*Y(:,i);
    end
end
multi = isfield(opts, 'Xtr');                 % network: lags of each bus net load
if multi
    F = price_features_multi(opts.Xtr, hrs, memH);
else
    F = price_features(xtr, hrs, memH);
    if ~useHarm, F = F(:, 1:memH); end
end
idx = (memH*4+1):n;                           % full history available
mu = mean(F(idx,:)); sd = std(F(idx,:)) + eps;
Z = [ones(numel(idx),1), (F(idx,:) - mu)./sd];
lam = 0.1*numel(idx);
B = (Z'*Z + lam*diag([0, ones(1, size(Z,2)-1)])) \ (Z'*Y(idx,:));
if isfield(opts, 'constant') && opts.constant
    B = [mean(Y(idx,:)); zeros(size(B,1)-1, N)];
end
pp = struct('B', B, 'mu', mu, 'sd', sd, 'memH', memH, 'useHarm', useHarm, 'cal', cal, 'multi', multi, ...
    'Omega', sqrt(pol.kappa./pol.ceff), 'meanPrice', mean(Ht.price), ...
    'r2', 1 - sum((Y(idx,:) - Z*B).^2)./sum((Y(idx,:) - mean(Y(idx,:))).^2));
end
