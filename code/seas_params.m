% SEAS_PARAMS  Fleet and design parameters of the seasonal 14-bus study, for
% the SpectraLMP README.
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - mean(meta.xb, 1); nTrain = 60*spd;
net = net14_setup(inf(20, 1));
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005;
par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
G = load(fullfile(root, 'results', 'pol_s.mat')); pol = G.pol;
fprintf('std(xs) all %.2f train %.2f\n', std(xs), std(xs(1:nTrain)));
fprintf('P   %s\nE   %s\netac %s\netad %s\nc0  %s\nc   %s\n', mat2str(par.P), mat2str(par.E), mat2str(par.etac), mat2str(par.etad), mat2str(par.c0), mat2str(par.c, 4));
fprintf('ceff %s\nkappa %s\ntau %s\nq %.4f eps %.2f e0/E %.2f\n', mat2str(pol.ceff(:)', 4), mat2str(pol.kappa(:)', 4), mat2str(sqrt(pol.ceff(:)./pol.kappa(:))', 4), par.q, par.eps_risk, par.e0(1)/par.E(1));
fprintf('genA %s\ngenQ %s\ngenBus %s\n', mat2str(net.genA(:)', 3), mat2str(net.genQ(:)', 4), mat2str(net.genBus(:)'));
