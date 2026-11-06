% FIG_IDEA_S  Seasonal version of the fig_idea cache: hindsight LMP at bus 1
% of the IEEE 14-bus case without a line limit, with the net load keeping
% its seasons (only its mean removed) and the energy drawn from the initial
% stock scored, as in netcfg sets ending in '14s'. Saves
% results/fig_idea_s.mat (lam, tau, dt) for fig_frame1c and fig_results6.
gp = getenv('GUROBI_HOME'); if ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end   % set GUROBI_HOME for other installs
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - mean(meta.xb, 1); nTrain = 60*spd;
net = net14_setup(inf(20, 1)); net.unitBus = [3 4 13 2]; net.genA(1:2) = 0.30;
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2); Xtest = X(nTrain+1:end, :);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005;
par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs; par.termval = true;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
ev = evalc('H = method_hindsight_net(Xtest, par, net, struct(''terminalFree'', false));'); %#ok<NASGU>
lam = H.lmp(:, 1); tau = sqrt(pol.ceff(:) ./ pol.kappa(:));
save(fullfile(root, 'results', 'fig_idea_s.mat'), 'lam', 'tau', 'dt');
fprintf('[fig_idea_s] tau = %s h, lb %.1f\n', mat2str(round(tau', 2)), H.cost);
