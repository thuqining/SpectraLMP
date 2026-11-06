% FIG_BODE14  Letter Fig. 1, right panel: filter bank (5) of the IEEE 14-bus
% fleet with crossovers (dotted) and bands (shaded), in Times New Roman 8 pt.
% Same design as net14_params.m, with each curve labeled at its peak in
% place of a legend, so the panel is shorter. Saves Figure/fig_bode.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xb = meta.xb - movmean(meta.xb, [7*spd-1, 0], 1); nTrain = 60*spd;
Fv = inf(20, 1); Fv(10) = 3.1;
net = net14_setup(Fv); net.unitBus = [3 4 13 2];
X = net14_loads(Xb, net, nTrain); xs = sum(X, 2);
par = params_fleet(xs, dt); par.a_fuel = 0.30; par.qs = 0.005; par.q = par.a_fuel/(2*std(xs(1:nTrain))) + par.qs;
ev = evalc('pol = synth_policy(xs(1:nTrain), par);'); %#ok<NASGU>
wc = zeros(1, par.N-1);
for i = 1:par.N-1
    wc(i) = sqrt((pol.kappa(i) - pol.kappa(i+1)) / (pol.ceff(i+1) - pol.ceff(i))) / (2*pi) * 24;
end
col = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;
gry = [0.45 0.45 0.45]; fn = 'Times New Roman'; fs = 8; durlab = {'1 h', '4 h', '30 h', '200 h'};
fcd = pol.f * 24;
fb = figure('Units', 'inches', 'Position', [0.5 0.5 2.67 1.28], 'Color', 'w', 'Visible', 'off');
axb = axes(fb); hold(axb, 'on'); box(axb, 'on');
edges = [max(fcd), wc, min(fcd)];
for i = 1:par.N
    patch(axb, [edges(i+1) edges(i) edges(i) edges(i+1)], [0 0 0.66 0.66], col(i,:), ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for i = 1:par.N, plot(axb, fcd, pol.H(i,:), '-', 'Color', col(i,:), 'LineWidth', 1.4); end
for i = 1:par.N-1, xline(axb, wc(i), ':', 'Color', gry, 'LineWidth', 0.8); end
% label each curve just above its peak
for i = 1:par.N
    [hm, k] = max(pol.H(i,:)); ha = 'center';
    if k <= 3, ha = 'left'; elseif k >= numel(fcd) - 2, ha = 'right'; end
    text(axb, fcd(k), hm + 0.06, durlab{i}, 'Color', col(i,:), 'FontName', fn, 'FontSize', fs, ...
        'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
end
set(axb, 'XScale', 'log', 'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, ...
    'XLim', [min(fcd) max(fcd)], 'YLim', [0 0.66], 'YTick', [0 0.2 0.4 0.6], 'XTick', [1e-2 1e-1 1 10], 'TickDir', 'out');
xlabel(axb, 'Frequency (cycles/day)'); ylabel(axb, 'Share of net load');
exportgraphics(fb, fullfile(root, 'Figure', 'fig_bode.pdf'), 'ContentType', 'vector'); close(fb);
fprintf('[fig_bode14] crossovers %s cycles/day\n', mat2str(wc, 3));
