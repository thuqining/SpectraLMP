% NET14_GAP  Predicted loss from causality (LQ model, finite-memory causal
% FIR vs two-sided FIR, pol.lqGap) on the IEEE 14-bus case without line
% limits, in points of the hindsight optimum of (1), for comparison with
% the realized causality loss (Alg. 1 with predicted vs exact averages).
gp = getenv('GUROBI_HOME'); if ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end   % set GUROBI_HOME for other installs
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - movmean(meta.xb, [7*spd-1, 0], 1); nTrain = 60*spd;
Fv = inf(20, 1);
net = net14_setup(Fv); net.unitBus = [3 4 13 2];
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2); Xtest = X(nTrain+1:end, :); T = size(Xtest, 1);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005; par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
opt = struct('terminalFree', true);
ev = evalc('H0 = method_hindsight_net(Xtest, par, net, opt);'); %#ok<NASGU>
gapUSD = pol.lqGap * dt * T;
fprintf('lqGap per h = %.4g, over test = %.1f $, hindsight (free) = %.1f $\n', pol.lqGap, gapUSD, H0.cost);
fprintf('Predicted causality loss = %.2f points of the hindsight optimum\n', 100*gapUSD/H0.cost);
fprintf('Jc_pred = %.4g, Jnc_pred = %.4g per h\n', pol.Jc_pred, pol.Jnc_pred);
save(fullfile(root, 'results', 'net14_gap.mat'), 'pol', 'gapUSD');
