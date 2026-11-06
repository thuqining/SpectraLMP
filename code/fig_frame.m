% FIG_FRAME  Letter Fig. 1 panels in schematic style (no ticks): one week
% of the uncongested IEEE 14-bus case. Offline panel: hindsight LMP and its
% two-sided averages (4). Online panels: recent net load (input) and SoC of
% the four units under Algorithm 1 (output; results/net/i14_frame_free.mat
% from fig_frame_run). Boxes, arrows, and labels are drawn in TikZ in the
% letter at the axes positions below (cm from the lower left of the page).
% Saves Figure/fig_idea.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
R = load(fullfile(root, 'results', 'net', 'i14_frame_free.mat')); tr = R.res.traj;
C = load(fullfile(root, 'results', 'fig_idea.mat'), 'lam', 'tau', 'dt'); dt = C.dt; spd = 96;
lam = C.lam; T = numel(lam); f = (0:T-1)'/T; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
N = numel(C.tau); mu = zeros(T, N);
for i = 1:N, mu(:, i) = real(ifft(fft(lam) ./ (1 + C.tau(i)^2 * w.^2))); end
win = 7*spd; best = 1; bv = -inf;                   % the week where the 30-h average varies most
for s = 30*spd:spd:T - 30*spd - win
    v = std(mu(s:s+win-1, 3)); if v > bv, bv = v; best = s; end
end
idx = best:best+win-1; td = (0:win-1)*dt/24;
x = sum(tr.Xtest(idx, :), 2); soc = tr.e(:, idx)' ./ tr.E(:)';
col = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;  % Okabe-Ito, as the filter bank
gry = [0.55 0.55 0.55]; lab = {'1 h', '4 h', '30 h', '200 h'};
fn = 'Times New Roman'; fs = 8; Wc = 11.7; Hc = 3.60; cm = 1/2.54;
AX = [4.55 2.10 2.75 1.04;                          % offline: LMP and averages
      0.70 0.15 2.00 1.12;                          % online: net load
      7.35 0.15 2.10 1.12];                         % online: SoC
fg = figure('Units', 'inches', 'Position', [0.5 0.5 Wc*cm Hc*cm], 'Color', 'none', 'InvertHardcopy', 'off', 'Visible', 'off');   % transparent: TikZ bands show
st = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1, ...
    'TitleFontWeight', 'normal', 'XLim', [0 7], 'XTick', [], 'YTick', [], 'Box', 'on', ...
    'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'Color', 'w'};
ax = gobjects(1, 3);
for k = 1:3
    ax(k) = axes(fg, 'Position', AX(k,:) ./ [Wc Hc Wc Hc]); hold(ax(k), 'on');
end
plot(ax(1), td, lam(idx), '-', 'Color', gry, 'LineWidth', 0.6);
for i = 1:N, plot(ax(1), td, mu(idx, i), '-', 'Color', col(i,:), 'LineWidth', 1.2); end
plot(ax(2), td, x, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.6);
for i = 1:N, plot(ax(3), td, soc(:, i), '-', 'Color', col(i,:), 'LineWidth', 1.0); end
yl = [min(lam(idx)) max(lam(idx))];
set(ax(1), st{:}, 'YLim', yl + [-0.05 0.05]*diff(yl));
set(ax(2), st{:}, 'YLim', [min(x) max(x)] + [-0.08 0.08]*range(x));
set(ax(3), st{:}, 'YLim', [-0.03 1.03]);
tt = {'LMP (gray) and its average (4) for each unit', 'Net load', 'SoC of each unit'};
for k = 1:3
    title(ax(k), tt{k}, 'FontSize', fs);                % time runs left to right
end
% direct labels at the right ends of the averages and SoC, spread apart
ya = mu(idx(end), :); [ys, o] = sort(ya); g = 0.24*diff(yl);
for k = 2:4, ys(k) = max(ys(k), ys(k-1) + g); end
ys = ys - max(0, ys(end) - yl(2));
for k = 1:4
    text(ax(1), 7.15, ys(k), lab{o(k)}, 'Color', col(o(k),:), 'FontName', fn, 'FontSize', fs - 1, ...
        'VerticalAlignment', 'middle', 'Clipping', 'off');
end
ye = soc(end, :); [ys, o] = sort(ye); for k = 2:4, ys(k) = max(ys(k), ys(k-1) + 0.27); end
ys = ys - max(0, ys(end) - 0.98);
for k = 1:4
    text(ax(3), 7.2, ys(k), lab{o(k)}, 'Color', col(o(k),:), 'FontName', fn, 'FontSize', fs, ...
        'VerticalAlignment', 'middle', 'Clipping', 'off');
end
% print keeps the page at Wc x Hc (no cropping), so the TikZ overlay can use AX
set(fg, 'PaperUnits', 'inches', 'PaperSize', [Wc Hc]*cm, 'PaperPosition', [0 0 Wc Hc]*cm);
print(fg, fullfile(root, 'Figure', 'fig_idea.pdf'), '-dpdf', '-vector');
close(fg);
fprintf('[fig_frame] week from day %.1f | SoC end %s\n', best*dt/24, mat2str(ye, 2));
