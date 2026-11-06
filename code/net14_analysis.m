% NET14_ANALYSIS  Network-specific results of the IEEE 14-bus case (line 5-6
% at 3.1 kW, base placement 1 h @ 3, 4 h @ 4, 30 h @ 13, 200 h @ 2):
%  (1) Theorem 1 on problem (1): R^2 of each unit's SoC dual explained by
%      the time-scale average (3) of the LMP at its bus, and of the system
%      (slack) price, with and without the line limit.
%  (2) Variance share by band of the system price and of the congestion
%      component (LMP at bus 13 minus slack price).
%  (3) Congestion and LMP spread without storage.
%  (4) Placement: oracle cost with the 30-h and 200-h units swapped and with
%      all units at bus 2.
% Saves results/net14_analysis.mat.
gp = getenv('GUROBI_HOME'); if ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end   % set GUROBI_HOME for other installs
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - movmean(meta.xb, [7*spd-1, 0], 1); nTrain = 60*spd;
Fv = inf(20, 1); Fv(10) = 3.1;
net = net14_setup(Fv); net.unitBus = [3 4 13 2]; net0 = net; net0.F(:) = inf;
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2); Xtest = X(nTrain+1:end, :); T = size(Xtest, 1);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005; par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
opt = struct('terminalFree', true);
ev = evalc('Ho = method_hindsight_net(Xtest, par, net, opt);'); %#ok<NASGU>
ev = evalc('H0 = method_hindsight_net(Xtest, par, net0, opt);'); %#ok<NASGU>

f = (0:T-1)'/T; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
avg = @(y, i) real(ifft(fft(y(:)) .* (pol.kappa(i)./(pol.kappa(i) + pol.ceff(i)*w.^2))));
r2 = @(D, y) 1 - sum((D - [ones(T,1), y]*([ones(T,1), y]\D)).^2)/sum((D - mean(D)).^2);
R2 = zeros(par.N, 3);                                   % [bus LMP, system price] limited; bus LMP free
for i = 1:par.N
    R2(i,1) = r2(Ho.shadow_value(i,:)', avg(Ho.priceUnit(:,i), i));
    R2(i,2) = r2(Ho.shadow_value(i,:)', avg(Ho.lmp(:,1), i));
    R2(i,3) = r2(H0.shadow_value(i,:)', avg(H0.priceUnit(:,i), i));
end

per = 2*pi./abs(w)/24;                                  % period, days
bands = {per > 7, per <= 7 & per > 1, per <= 1 & per > 4/24, per <= 4/24 & isfinite(per)};
P = @(y) abs(fft(y(:) - mean(y))).^2;
sh = @(y) cellfun(@(b) sum(P(y).*b)/sum(P(y)), bands);
cong = Ho.lmp(:,13) - Ho.lmp(:,1);
A.shSys = sh(Ho.lmp(:,1)); A.shCong = sh(cong); A.stdSys = std(Ho.lmp(:,1)); A.stdCong = std(cong);

parNS = par; parNS.P = 1e-6*par.P; parNS.E = 1e-6*par.E; parNS.e0 = 1e-6*par.e0; parNS.eref = parNS.e0;
ev = evalc('Hn = method_hindsight_net(Xtest, parNS, net, opt);'); %#ok<NASGU>
spread = @(H) mean(max(H.lmp,[],2) - min(H.lmp,[],2));
pl = {[3 4 2 13], [2 2 2 2]}; Cp = zeros(1, 2); Gp = Cp;
for k = 1:2
    nk = net; nk.unitBus = pl{k};
    ev = evalc('Hk = method_hindsight_net(Xtest, par, nk, opt);'); %#ok<NASGU>
    Cp(k) = Hk.cost; Gp(k) = Hk.congested;
end
A.R2 = R2; A.oracle = Ho.cost; A.oracleFree = H0.cost; A.cong = Ho.congested; A.spread = spread(Ho);
A.congNS = Hn.congested; A.spreadNS = spread(Hn); A.costNS = Hn.cost; A.place = Cp; A.placeCong = Gp;
save(fullfile(root, 'results', 'net14_analysis.mat'), 'A');
fprintf('(1) R^2 of SoC dual: [bus LMP avg | system price avg] limited, [bus LMP avg] free\n');
for i = 1:par.N
    fprintf('   unit %d (%gh @ bus %d): %.3f | %.3f , free %.3f\n', i, par.h(i), net.unitBus(i), R2(i,:));
end
fprintf('(2) band shares [>1wk 1d-1wk 4h-1d <4h]: system %s (std %.3f) | congestion %s (std %.3f)\n', ...
    mat2str(round(A.shSys,2)), A.stdSys, mat2str(round(A.shCong,2)), A.stdCong);
fprintf('(3) congestion %.1f%% with storage, %.1f%% without; spread %.3f vs %.3f; cost %.0f vs %.0f $ (free %.0f)\n', ...
    100*A.cong, 100*A.congNS, A.spread, A.spreadNS, A.oracle, A.costNS, A.oracleFree);
fprintf('(4) swap 30h<->200h: %+.1f%% (congested %.1f%%); all at bus 2: %+.1f%% (congested %.1f%%)\n', ...
    100*(Cp(1)/Ho.cost-1), 100*Gp(1), 100*(Cp(2)/Ho.cost-1), 100*Gp(2));
