% SEAS_NUMBERS  Case-study numbers of the seasonal 14-bus study (sets '14s'):
%  (1) Theorem 1: R^2 of each unit's hindsight SoC dual explained by the
%      average (4) of the LMP at its bus, without and with the line limit;
%  (2) Corollary 1: empirical frequency split of the proposed dispatch
%      (i14s_main_free) against the shares H_i of (5), saved to
%      results/spectral_share_s.mat for fig_results6.
gp = getenv('GUROBI_HOME'); if ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end   % set GUROBI_HOME for other installs
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - mean(meta.xb, 1); nTrain = 60*spd;
Fv = inf(20, 1); Fv(10) = 3.3;
net = net14_setup(Fv); net.unitBus = [3 4 13 2]; net.genA(1:2) = 0.30; net0 = net; net0.F(:) = inf;
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2); Xtest = X(nTrain+1:end, :); T = size(Xtest, 1);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005;
par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs; par.termval = true;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
opt = struct('terminalFree', false);
ev = evalc('Ho = method_hindsight_net(Xtest, par, net, opt);'); %#ok<NASGU>
ev = evalc('H0 = method_hindsight_net(Xtest, par, net0, opt);'); %#ok<NASGU>
f = (0:T-1)'/T; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
avg = @(y, i) real(ifft(fft(y(:)) .* (pol.kappa(i)./(pol.kappa(i) + pol.ceff(i)*w.^2))));
r2 = @(D, y) 1 - sum((D - [ones(T,1), y]*([ones(T,1), y]\D)).^2)/sum((D - mean(D)).^2);
R2 = zeros(par.N, 2);
for i = 1:par.N
    R2(i,1) = r2(H0.shadow_value(i,:)', avg(H0.priceUnit(:,i), i));
    R2(i,2) = r2(Ho.shadow_value(i,:)', avg(Ho.priceUnit(:,i), i));
end
fprintf('tau = %s h\n', mat2str(round(sqrt(pol.ceff(:)./pol.kappa(:))', 2)));
fprintf('lb free %.1f lim %.1f congested %.1f%%\n', H0.cost, Ho.cost, 100*Ho.congested);
fprintf('R2 [free lim] per unit:\n'); disp(round(R2, 3));
fprintf('hindsight mu std/mean per unit (free): %s\n', mat2str(round(std(H0.shadow_value, 0, 2)'./mean(H0.shadow_value, 2)', 3)));
% (2) frequency split
S = load(fullfile(root, 'results', 'net', 'i14s_splitmain_free.mat')); tr = S.res.traj;
x = sum(tr.Xtest, 2); P = tr.p; if size(P,1) ~= numel(x), P = P'; end
c = pol.ceff(:); k = pol.kappa(:); q = par.q;
nw = 96*30; win = hann(nw); nov = nw/2;
[Sxx, fr] = cpsd(x - mean(x), x - mean(x), win, nov, nw, 1/dt);
cyc = fr*24; He = zeros(numel(fr), size(P, 2));
for i = 1:size(P, 2)
    Spx = cpsd(P(:,i) - mean(P(:,i)), x - mean(x), win, nov, nw, 1/dt);
    He(:,i) = real(Spx ./ Sxx);
end
ww = 2*pi*fr; a = c' + k'./max(ww,1e-9).^2; Ht = (1./a) ./ (1/q + sum(1./a, 2));
sel = cyc >= 0.05 & cyc <= 40;
for b = [0.1 0.3 1 3 10 30]
    [~, j] = min(abs(cyc - b));
    fprintf('%5.1f c/d | emp %s | thy %s\n', cyc(j), mat2str(round(He(j,:),2)), mat2str(round(Ht(j,:),2)));
end
save(fullfile(root, 'results', 'spectral_share_s.mat'), 'cyc', 'He', 'Ht', 'sel');
fprintf('q %.4f\n', q);
