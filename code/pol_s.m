% POL_S  Calibrated weights (Step 1) of the seasonal base fleet, saved to
% results/pol_s.mat for the theory curves of Fig. 2(c).
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - mean(meta.xb, 1); nTrain = 60*spd;
net = net14_setup(inf(20, 1));
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005;
par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
save(fullfile(root, 'results', 'pol_s.mat'), 'pol');
fprintf('tau = %s h\n', mat2str(round(sqrt(pol.ceff(:)./pol.kappa(:))', 2)));
