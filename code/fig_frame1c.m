% FIG_FRAME1C  Data panels of Fig. 1 (single column, 8.85 x 3.9 cm): one
% week of the uncongested IEEE 14-bus case, with time axes in days.
% Offline: hindsight LMP (gray) and its averages (4) for the four units.
% Online: SoC of the four units under Algorithm 1
% (results/net/i14_frame_free.mat from fig_frame_run). Boxes, arrows, and
% titles are added by Figure/fig_framework.tex at the axes positions AX
% below (cm from the lower left). Saves Figure/fig_frame1c.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
R = load(fullfile(root, 'results', 'net', 'i14s_splitmain_free.mat')); tr = R.res.traj;
C = load(fullfile(root, 'results', 'fig_idea_s.mat'), 'lam', 'tau', 'dt'); dt = C.dt; spd = 96;
lam = 1000*C.lam; T = numel(lam);                   % LMP in $/MWh
f = (0:T-1)'/T; f(f > 0.5) = f(f > 0.5) - 1; w = 2*pi*f/dt;
N = numel(C.tau); mu = zeros(T, N);
for i = 1:N, mu(:, i) = real(ifft(fft(lam) ./ (1 + C.tau(i)^2 * w.^2))); end
if ~exist('NDAYS', 'var'), NDAYS = 28; end
win = NDAYS*spd; best = 1; bv = -inf;               % the window where the 30-h average varies most
for s = 30*spd:spd:T - 30*spd - win
    v = std(mu(s:s+win-1, 3)); if v > bv, bv = v; best = s; end
end
if ~exist('S0DAY', 'var'), S0DAY = 17; end                % autumn week: the four SoC patterns differ most (pick_week.m)
best = round(S0DAY*spd) + 1;
idx = best:best+win-1; td = (0:win-1)*dt/24;
soc = tr.e(:, idx)' ./ tr.E(:)';
col = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;
gry = [0.55 0.55 0.55]; lab = {'1 h', '4 h', '30 h', '200 h'};
fn = 'Times New Roman'; fs = 6.5; Wc = 8.85; Hc = 3.9; cm = 1/2.54;
AX = [3.85 2.58 4.22 0.94;                          % offline: LMP and averages
      3.85 0.58 4.22 0.94];                         % online: SoC
fg = figure('Units', 'inches', 'Position', [0.5 0.5 Wc*cm Hc*cm], 'Color', 'none', ...
    'InvertHardcopy', 'off', 'Visible', 'off');
st = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'XLim', [0 NDAYS], 'XTick', 0:max(1, NDAYS/7):NDAYS, ...
    'Box', 'on', 'TickDir', 'in', 'TickLength', [0.012 0.012], 'XColor', [0.25 0.25 0.25], ...
    'YColor', [0.25 0.25 0.25], 'Color', 'w', 'Layer', 'top'};
ax1 = axes(fg, 'Position', AX(1,:) ./ [Wc Hc Wc Hc]); hold(ax1, 'on');
ax2 = axes(fg, 'Position', AX(2,:) ./ [Wc Hc Wc Hc]); hold(ax2, 'on');
plot(ax1, td, lam(idx), '-', 'Color', gry, 'LineWidth', 0.6);
for i = 1:N, plot(ax1, td, mu(idx, i), '-', 'Color', col(i,:), 'LineWidth', 1.0); end
for i = 1:N, plot(ax2, td, soc(:, i), '-', 'Color', col(i,:), 'LineWidth', 0.9); end
yl = [min(lam(idx)) max(lam(idx))]; yt = 0:200:yl(2);
set(ax1, st{:}, 'YLim', yl + [-0.05 0.05]*diff(yl), 'YTick', yt);
set(ax2, st{:}, 'YLim', [-0.03 1.03], 'YTick', [0 0.5 1]);
xlabel(ax1, 'Time (day)'); xlabel(ax2, 'Time (day)');
for a = [ax1 ax2], a.XLabel.Units = 'normalized'; a.XLabel.Position(2) = -0.33; end
% direct labels right of each panel, spread apart; LMP in gray
ya = [mu(idx(end), :) lam(idx(end))]; lb = [lab {'LMP'}]; cl = [col; gry*0.8];
[ys, o] = sort(ya); g = 0.225*diff(yl);
for k = 2:5, ys(k) = max(ys(k), ys(k-1) + g); end
ys = ys - max(0, ys(end) - yl(2) - 0.04*diff(yl));
for k = 1:5
    text(ax1, NDAYS*1.014, ys(k), lb{o(k)}, 'Color', cl(o(k),:), 'FontName', fn, 'FontSize', fs, ...
        'VerticalAlignment', 'middle', 'Clipping', 'off');
end
ye = soc(end, :); [ys, o] = sort(ye); ys(1) = max(ys(1), 0.12);   % keep clear of the x tick label
for k = 2:4, ys(k) = max(ys(k), ys(k-1) + 0.26); end
ys = ys - max(0, ys(end) - 1.0);
for k = 1:4
    text(ax2, NDAYS*1.014, ys(k), lab{o(k)}, 'Color', col(o(k),:), 'FontName', fn, 'FontSize', fs, ...
        'VerticalAlignment', 'middle', 'Clipping', 'off');
end
set(fg, 'PaperUnits', 'inches', 'PaperSize', [Wc Hc]*cm, 'PaperPosition', [0 0 Wc Hc]*cm);
if ~exist('OUTNAME', 'var'), OUTNAME = 'fig_frame1c.pdf'; end
print(fg, fullfile(root, 'Figure', OUTNAME), '-dpdf', '-vector');
close(fg);
fprintf('[fig_frame1c] week from day %.1f, LMP ticks %s\n', best*dt/24, mat2str(yt));
