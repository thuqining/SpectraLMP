function fleet_refit(histDays, cfgs, recal, calMin, minHist)
% FLEET_REFIT  Fleet study under the deployment protocol: the proposed
% policy and ADP start with histDays of history before the test year and
% are refitted every 30 days on all data since then. MPC and the myopic
% rule need no training and are run for histDays = 60 only.
% v27 accounting: the oracle is problem (1) without any terminal condition
% (terminalFree), and every method is scored on the operating cost (1a)
% over the test year (res.prop etc.); the settled costs are kept as
% res.*Settled. With the default options (the adopted algorithm, v26c:
% constant values until 16 days of data) it saves
% results/fleet_v27_<histDays>_<k>.mat.
% Experimental options (v26, rejected): recal = monthly kappa
% recalibration; calMin = duration prior until calMin days. They save
% fleet_recal_* / fleet_cal<calMin>_* instead.
if nargin < 3, recal = false; end
if nargin < 4, calMin = 0; end
if nargin < 5, minHist = 16; end
gp = getenv('GUROBI_HOME'); if ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end   % set GUROBI_HOME for other installs
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[xall, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
hrs = mod(hour(meta.t0) + minute(meta.t0)/60 + (0:numel(xall)-1)'*dt, 24);
xs = xall - movmean(xall, [7*spd-1, 0]);
nTrain = 60*spd; xtest = xs(nTrain+1:end);
i0 = nTrain - histDays*spd + 1;                     % commissioning
xsub = xs(i0:end); hsub = hrs(i0:end); n0 = histDays*spd;
Ptot = 137.8;
CFG = { 'N=1', [8]; 'N=2 close', [4, 30]; 'N=2 wide', [1, 200];
        'N=4 close', [4, 7.8, 15.4, 30]; 'N=4 wide', [1, 4, 30, 200];
        'N=6 close', [4, 6, 9, 13.4, 20, 30]; 'N=6 wide', [1, 2.9, 8.4, 24, 69, 200] };
hA = [1, 4, 30, 200]; c0A = [0.001, 0.003, 0.010, 0.030];
etcA = [0.96, 0.95, 0.85, 0.65]; etdA = [0.96, 0.95, 0.85, 0.55];
cfun = @(h) exp(interp1(log(hA), log(c0A), log(min(max(h,1),200)), 'linear'));
etcf = @(h) interp1(log(hA), etcA, log(min(max(h,1),200)), 'linear');
etdf = @(h) interp1(log(hA), etdA, log(min(max(h,1),200)), 'linear');
adopted = ~recal && calMin == 0 && minHist == 16;
opc = @(R) R.cost + R.settlement;                   % operating cost (1a)
for k = cfgs
    h = CFG{k,2}; N = numel(h); w = 0.75.^(0:N-1); w = w/sum(w);
    par = params_fleet(xs, dt);
    par.qs = 0.005; par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
    par.N = N; par.h = h; par.P = Ptot*w; par.E = par.P.*h;
    par.c0 = cfun(h); par.etac = etcf(h); par.etad = etdf(h);
    par.c = par.c0 + 0.01*(1 - par.etac.*par.etad);
    par.gamma = 0.1./h; par.eref = 0.6*par.E; par.e0 = par.eref;
    par.name = arrayfun(@(x) sprintf('%gh', x), h, 'UniformOutput', false);
    ev = evalc('pol = synth_policy(xsub(1:n0), par);'); %#ok<NASGU>
    vspec = struct('pol', pol, 'opts', struct('lossCal', true), 'recal', recal);
    if calMin > 0, vspec = struct('pol', pol, 'opts', struct('lossCal', true), 'calMin', calMin); end
    if minHist > 0, vspec.minHist = minHist; end
    out = refit_values(xsub, hsub, par, n0, 30, {vspec});
    ev = evalc('Rp = method_price(xtest, xsub(1:n0), hsub, par, [], ''ext'', out.V{1});'); %#ok<NASGU>
    Ra = run_adp_blocks(xtest, par, out.fits, out.blk);
    ev = evalc('Rh = method_hindsight(xtest, setfield(par, ''terminalFree'', true));'); %#ok<NASGU>
    res = struct('cfg', CFG{k,1}, 'lb', Rh.cost, 'prop', opc(Rp), 'adp', opc(Ra), 'mpc', NaN, 'myop', NaN, ...
        'propSettled', Rp.cost, 'adpSettled', Ra.cost, 'acct', 'v27 operating cost, free oracle');
    if histDays == 60 && adopted
        ev = evalc('Rm = method_mpc(xtest, par, pol, ''persistence'');'); %#ok<NASGU>
        ev = evalc('Ry = method_myopic(xtest, par, pol);'); %#ok<NASGU>
        res.mpc = opc(Rm); res.myop = opc(Ry);
    end
    res.tau = cellfun(@(q) sqrt(q.ceff./q.kappa), out.pol, 'UniformOutput', false);
    if adopted, tag = 'fleet_v27';
    elseif recal, tag = 'fleet_recal';
    elseif calMin > 0, tag = sprintf('fleet_cal%d', calMin);
    else, tag = sprintf('fleet_mh%d', minHist);
    end
    save(fullfile(root, 'results', sprintf('%s_%d_%d.mat', tag, histDays, k)), 'res');
    fprintf('[fleet %d d] %-10s Prop %.1f ADP %.1f MPC %.1f Myop %.1f\n', histDays, CFG{k,1}, ...
        100*res.prop/res.lb, 100*res.adp/res.lb, 100*res.mpc/res.lb, 100*res.myop/res.lb);
end
end
