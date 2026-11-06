% FIG_DECOMP  Letter Fig. 1 (right): the hindsight LMP of the uncongested
% IEEE 14-bus case against the two-sided average (4) of each unit, over a
% window matched to its time constant (2 days, 1 week, 2 months, 1 year,
% nested around the week of fig_frame.m). Red: LMP above the average
% (discharge); blue: below (charge). Long windows show the LMP in bins
% (15 min, 15 min, 2 h, 1 day). The title is drawn in TikZ.
% Saves Figure/fig_decomp.pdf (6.1 x 3.6 cm).
here = fileparts(mfilename('fullpath')); root = fileparts(here);
C = load(fullfile(root, 'results', 'fig_idea.mat'), 'lam', 'tau', 'dt'); dt = C.dt; spd = 96;
lam = C.lam(:); T = numel(lam); f = (0:T-1)'/T; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
N = numel(C.tau); mu = zeros(T, N);
for i = 1:N, mu(:, i) = real(ifft(fft(lam) ./ (1 + C.tau(i)^2 * w.^2))); end
win = 7*spd; best = 1; bv = -inf;                   % same week as fig_frame.m
for s = 30*spd:spd:T - 30*spd - win
    v = std(mu(s:s+win-1, 3)); if v > bv, bv = v; best = s; end
end
mid = best + win/2;
len = [2 7 60 365]*spd; bin = [1 1 8 96];           % window and bin (steps) per unit
col = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;
cD = [214 96 77]/255; cC = [67 147 195]/255;        % discharge, charge
lab = {'1 h', '4 h', '30 h', '200 h'}; wl = {'2 days', '1 week', '2 months', '1 year'};
fn = 'Times New Roman'; fs = 8; Wc = 6.1; Hc = 3.6; cm = 1/2.54;
x0 = 1.15; wa = 4.85; ha = 0.62; gap = 0.10; y0 = 0.08;
fg = figure('Units', 'inches', 'Position', [0.5 0.5 Wc*cm Hc*cm], 'Color', 'w', 'Visible', 'off');
for i = 1:N
    s = max(1, min(T - len(i) + 1, round(mid - len(i)/2))); idx = s:s+len(i)-1;
    nb = floor(len(i)/bin(i)); idx = idx(1:nb*bin(i));
    L = mean(reshape(lam(idx), bin(i), nb), 1)'; M = mean(reshape(mu(idx, i), bin(i), nb), 1)';
    t = ((1:nb)' - 0.5)/nb;
    ax = axes(fg, 'Position', [x0, y0 + (N - i)*(ha + gap), wa, ha] ./ [Wc Hc Wc Hc]); hold(ax, 'on');
    up = max(L, M); dn = min(L, M);
    fill(ax, [t; flipud(t)], [up; flipud(M)], cD, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    fill(ax, [t; flipud(t)], [M; flipud(dn)], cC, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    plot(ax, t, L, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.4);
    plot(ax, t, M, '-', 'Color', col(i,:), 'LineWidth', 0.9);
    yl = [min(L) max(L)];
    set(ax, 'XLim', [0 1], 'YLim', yl + [-0.06 0.06]*diff(yl), 'XTick', [], 'YTick', [], 'Box', 'on', ...
        'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], 'Color', 'w');
    text(ax, -0.03, 0.62, lab{i}, 'Units', 'normalized', 'Color', col(i,:), 'FontName', fn, ...
        'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(ax, -0.03, 0.22, wl{i}, 'Units', 'normalized', 'Color', [0.35 0.35 0.35], 'FontName', fn, ...
        'FontSize', fs - 1, 'HorizontalAlignment', 'right');
end
set(fg, 'PaperUnits', 'inches', 'PaperSize', [Wc Hc]*cm, 'PaperPosition', [0 0 Wc Hc]*cm);
print(fg, fullfile(root, 'Figure', 'fig_decomp.pdf'), '-dpdf', '-vector');
close(fg);
fprintf('[fig_decomp] center day %.1f\n', mid*dt/24);
