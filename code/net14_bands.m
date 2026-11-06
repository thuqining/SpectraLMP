% NET14_BANDS  Crossovers of the cost spectra (Cor. 1) and the value of
% each frequency band (Cor. A1) for the IEEE 14-bus fleet without line
% limits, from the calibration saved by net14_gap.m.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
S = load(fullfile(root, 'results', 'net14_gap.mat')); pol = S.pol;
c = pol.ceff(:); k = pol.kappa(:); N = numel(c);
spd = 96; dt = 0.25;
xs = [];  %#ok<NASGU>
q = 0.0073;
% crossovers between neighbours, in cycles/day (omega in rad/h)
fprintf('tau (h): %s\n', mat2str(sqrt(c./k)', 3));
for i = 1:N-1
    w = sqrt((k(i)-k(i+1))/(c(i+1)-c(i)));
    fprintf('crossover %d-%d: %.3g cycles/day (period %.3g h)\n', i, i+1, w*24/(2*pi), 2*pi/w);
end
% band values: v = q(1-S)Phi_x on the FFT grid of the training autocovariance
ac = pol.ac(:); n = numel(ac); Phi = real(fft([ac; flipud(ac(2:end))]));
m = numel(Phi); f = (0:m-1)'/m; f(f > 0.5) = f(f > 0.5) - 1; w = abs(2*pi*f/dt);
w(w == 0) = eps;
a = c + k./(w'.^2);
Sw = (1/q) ./ (1/q + sum(1./a, 1))';
v = q*(1 - Sw).*max(Phi, 0);
per = 2*pi./w/24;
bands = {per > 7, per <= 7 & per > 1, per <= 1 & per > 4/24, per <= 4/24};
names = {'>1 wk', '1 d-1 wk', '4 h-1 d', '<4 h'};
tot = sum(v);
for b = 1:4
    fprintf('band %-9s share of storage savings %.1f%%\n', names{b}, 100*sum(v(bands{b}))/tot);
end
fprintf('storage saves %.1f%% of the generator-only LQ cost\n', 100*tot/sum(q*max(Phi,0)));
