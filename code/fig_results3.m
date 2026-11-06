% FIG_RESULTS3  Fig. 2, IEEE 14-bus system, in three rows (3.5 x 2.66 in):
% (a) hindsight LMP without congestion (gray) and the average (4) of each
%     unit (color) over a window matched to its time constant (2 days,
%     1 week, 2 months, 1 year, nested around the week of Fig. 1); red
%     where the LMP is above (discharge), blue below (charge);
% (b) optimality gap: extra cost over the hindsight optimum of the
%     proposed method split into the LQ approximation with the limits
%     (oracle) and the prediction error (Theorem 2), against ADP, MPC, and
%     MPC with perfect forecasts; bars without congestion, diamonds with it;
% (c) storage share of each unit by frequency: theory H_i/sum H (lines)
%     against the dispatch of Algorithm 1 over the test year (dots);
% (d) six fleets with 60 days of history: proposed and ADP joined;
%     circles without congestion, diamonds with it; MPC in text.
% Writes Figure/fig_results3.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
ld = @(t) getfield(load(fullfile(root, 'results', 'net', ['i14_' t '.mat'])), 'res'); %#ok<GFLD>
pc = @(r, f) 100*r.(f)/r.lb - 100;
fn = 'Times New Roman'; fs = 8;
toRGB = @(h) sscanf(h(2:end), '%2x%2x%2x')' / 255;
cP = toRGB('#B2182B'); cA = toRGB('#4D4D4D'); cM = toRGB('#B3B3B3');
cu = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;      % unit colors of Fig. 1
W = 3.5; H = 2.66; in2n = @(l, b, w, h) [l/W, b/H, w/W, h/H];
f = figure('Units', 'inches', 'Position', [0.5 0.5 W H], 'Color', 'w', 'Visible', 'off');
opt = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'TickDir', 'out', 'Box', 'off', ...
    'TickLength', [0.02 0.02]};

% ---- (a) LMP and the average of each unit ----
D = load(fullfile(root, 'results', 'fig_idea.mat'), 'lam', 'tau', 'dt'); spd = 96;
lam = D.lam(:); T = numel(lam); fq = (0:T-1)'/T; fq(fq > 0.5) = fq(fq > 0.5) - 1; wq = 2*pi*fq/D.dt;
mu = zeros(T, 4);
for i = 1:4, mu(:, i) = real(ifft(fft(lam) ./ (1 + D.tau(i)^2 * wq.^2))); end
win = 7*spd; best = 1; bv = -inf;                   % same week as Fig. 1
for s0 = 30*spd:spd:T - 30*spd - win
    v = std(mu(s0:s0+win-1, 3)); if v > bv, bv = v; best = s0; end
end
mid = best + win/2; len = [2 7 60 365]*spd; bin = [1 1 8 96];
cD = [214 96 77]/255; cC = [67 147 195]/255;
ul = {'1 h', '4 h', '30 h', '200 h'}; wl = {'2 d', '1 wk', '2 mo', '1 yr'};
ha = 0.150; ga = 0.025; ya0 = 1.79;
for i = 1:4
    s0 = max(1, min(T - len(i) + 1, round(mid - len(i)/2))); id = s0:s0+len(i)-1;
    nb = floor(len(i)/bin(i)); id = id(1:nb*bin(i));
    L = mean(reshape(lam(id), bin(i), nb), 1)'; M = mean(reshape(mu(id, i), bin(i), nb), 1)';
    tt = ((1:nb)' - 0.5)/nb;
    axa = axes(f, 'Position', in2n(0.60, ya0 + (4 - i)*(ha + ga), 2.82, ha)); hold(axa, 'on');
    up = max(L, M); dn = min(L, M);
    fill(axa, [tt; flipud(tt)], [up; flipud(M)], cD, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    fill(axa, [tt; flipud(tt)], [M; flipud(dn)], cC, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    plot(axa, tt, L, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.4);
    plot(axa, tt, M, '-', 'Color', cu(i,:), 'LineWidth', 0.9);
    yr = [min(L) max(L)];
    set(axa, 'XLim', [0 1], 'YLim', yr + [-0.06 0.06]*diff(yr), 'XTick', [], 'YTick', [], 'Box', 'on', ...
        'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], 'Color', 'w');
    text(axa, -0.115, 0.5, ul{i}, 'Units', 'normalized', 'Color', cu(i,:), 'FontName', fn, ...
        'FontSize', fs - 1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    text(axa, -0.012, 0.5, wl{i}, 'Units', 'normalized', 'Color', [0.35 0.35 0.35], 'FontName', fn, ...
        'FontSize', fs - 1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

% ---- (b) optimality gap ----
Mf = ld('main_free'); Ml = ld('main_lim'); Of = ld('hindcal_free'); Ol = ld('hindcal_lim');
Qf = ld('mpcfit_free'); Ql = ld('mpcfit_lim');
P = [pc(Mf,'prop'), pc(Ml,'prop')]; O = [pc(Of,'hindcal'), pc(Ol,'hindcal')];
A = [pc(Mf,'adp'), pc(Ml,'adp')]; Q = [pc(Qf,'mpcfit'), pc(Ql,'mpcfit')];
Qs = [pc(Mf,'mpcexact'), pc(Ml,'mpcexact')];
ax = axes(f, 'Position', in2n(0.36, 1.14, 1.30, 0.44)); hold(ax, 'on');
bw = 0.66; cL = cP*0.5 + 0.5;
bar(ax, 1, P(1), bw, 'FaceColor', cP, 'EdgeColor', 'none');
bar(ax, 1, O(1), bw, 'FaceColor', cL, 'EdgeColor', 'none');
bar(ax, 2, A(1), bw, 'FaceColor', cA, 'EdgeColor', 'none');
bar(ax, 3, Q(1), bw, 'FaceColor', cM, 'EdgeColor', 'none');
bar(ax, 4, Qs(1), bw, 'FaceColor', 'w', 'EdgeColor', cM, 'LineWidth', 0.8);
plot(ax, 1:4, [P(2) A(2) Q(2) Qs(2)], 'd', 'MarkerSize', 3, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'LineStyle', 'none');
text(ax, 1, O(1)/2, 'LQ', 'FontName', fn, 'FontSize', fs-2, 'Color', 'w', 'HorizontalAlignment', 'center');
text(ax, 1, (O(1)+P(1))/2 - 0.6, 'pred.', 'FontName', fn, 'FontSize', fs-2, 'Color', 'w', ...
    'HorizontalAlignment', 'center');
set(ax, opt{:}, 'XLim', [0.4 4.6], 'XTick', 1:4, 'XTickLabel', {'Prop.', 'ADP', 'MPC', 'MPC*'}, ...
    'XTickLabelRotation', 0, 'YLim', [0 22], 'YTick', 0:10:20, 'YGrid', 'on', ...
    'GridColor', [0.88 0.88 0.88], 'GridAlpha', 1);
ax.XAxis.FontSize = fs - 1;
ylabel(ax, 'Extra cost (%)');
text(ax, 4.55, 21.8, '(b)', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% ---- (c) storage share by frequency ----
S = load(fullfile(root, 'results', 'spectral_share.mat'));
cyc = S.cyc; sel = S.sel;
Ht = S.Ht ./ sum(S.Ht, 2); He = max(S.He, 0); He = He ./ sum(He, 2);
ax2 = axes(f, 'Position', in2n(2.12, 1.14, 1.30, 0.44)); hold(ax2, 'on');
edges = logspace(log10(0.05), log10(40), 14); cb = sqrt(edges(1:end-1).*edges(2:end));
for i = 1:4
    semilogx(ax2, cyc(sel), Ht(sel, i), '-', 'Color', cu(i,:), 'LineWidth', 1.0);
    hb = nan(size(cb));
    for k = 1:numel(cb)
        m = cyc >= edges(k) & cyc < edges(k+1); if any(m), hb(k) = mean(He(m, i)); end
    end
    semilogx(ax2, cb, hb, 'o', 'MarkerSize', 2.4, 'MarkerFaceColor', cu(i,:), 'MarkerEdgeColor', 'none');
end
set(ax2, opt{:}, 'XScale', 'log', 'XLim', [0.05 40], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'}, ...
    'YLim', [0 1], 'YTick', [0 0.5 1]);
ylabel(ax2, 'Storage share'); xlabel(ax2, 'Cycles/day');
text(ax2, 38, 1.0, '(c)', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
labs = {'1 h', '4 h', '30 h', '200 h'}; lx = [5.5 1.4 0.12 0.45]; ly = [0.84 0.72 0.82 0.95];
for i = 1:4
    text(ax2, lx(i), ly(i), labs{i}, 'Color', cu(i,:), 'FontName', fn, 'FontSize', fs-1);
end

% ---- (d) six fleets, 60 days ----
rel = fullfile(root, 'SpectraLMP', 'results'); if ~isfolder(rel), rel = fullfile(root, 'results'); end
T = readtable(fullfile(rel, 'paired_36.csv'), 'TextType', 'string');
T = T(T.history_days == 60, :);
fl = ["2c" "2w" "4c" "4w" "6c" "6w"]; lns = ["none" "3.1 kW"]; cap = 21;
ax3 = axes(f, 'Position', in2n(0.36, 0.32, 2.78, 0.42)); hold(ax3, 'on');
for k = 1:6
    mm = nan(1, 2);
    for l = 1:2
        r = T(T.fleet == fl(k) & T.line_limit == lns(l), :);
        y = k + (l - 1.5)*0.36; p = r.proposed - 100; a = r.adp - 100; mm(l) = r.mpc - 100;
        plot(ax3, [p a], [y y], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.2);
        mk = 'o'; if l == 2, mk = 'd'; end
        plot(ax3, a, y, mk, 'MarkerSize', 3.2, 'MarkerFaceColor', cA, 'MarkerEdgeColor', cA);
        plot(ax3, p, y, mk, 'MarkerSize', 3.2, 'MarkerFaceColor', cP, 'MarkerEdgeColor', cP);
    end
    text(ax3, cap + 0.2, k, sprintf('%.0f / %.0f', mm), 'FontName', fn, 'FontSize', fs-1, ...
        'VerticalAlignment', 'middle', 'Color', [0.35 0.35 0.35]);
end
text(ax3, cap + 0.2, 0.05, 'MPC', 'FontName', fn, 'FontSize', fs-1, 'VerticalAlignment', 'middle', ...
    'Color', [0.35 0.35 0.35]);
set(ax3, opt{:}, 'YLim', [0.4 6.6], 'YTick', 1:6, 'YTickLabel', cellstr(fl), 'YDir', 'reverse', ...
    'XLim', [4 cap], 'XTick', 5:5:20, 'XGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 1, ...
    'Clipping', 'off');
xlabel(ax3, 'Extra cost (%)'); ylabel(ax3, 'Fleet'); ax3.YAxis.FontSize = fs - 1;
% legend for (b) and (d), above both rows
hP = plot(ax3, nan, nan, 's', 'MarkerFaceColor', cP, 'MarkerEdgeColor', cP, 'MarkerSize', 3.6);
hA = plot(ax3, nan, nan, 's', 'MarkerFaceColor', cA, 'MarkerEdgeColor', cA, 'MarkerSize', 3.6);
hM = plot(ax3, nan, nan, 's', 'MarkerFaceColor', cM, 'MarkerEdgeColor', cM, 'MarkerSize', 3.6);
hD = plot(ax3, nan, nan, 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerSize', 3.2);
lg = legend(ax3, [hP hA hM hD], {'Proposed', 'ADP', 'MPC', '3.1-kW line limit'}, 'FontName', fn, ...
    'FontSize', fs, 'Box', 'off', 'Orientation', 'horizontal', 'NumColumns', 4);
lg.ItemTokenSize = [8 6]; lg.Position = in2n(0.42, 2.50, 3.06, 0.15);
% panel letters outside the axes
axo = axes(f, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 W], 'YLim', [0 H]); hold(axo, 'on');
pl = {'(a)', 0.02, 2.575; '(d)', 0.02, 0.82};
for k = 1:2
    text(axo, pl{k,2}, pl{k,3}, pl{k,1}, 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
end
exportgraphics(f, fullfile(root, 'Figure', 'fig_results3.pdf'), 'ContentType', 'vector'); close(f);
fprintf('gap: oracle %.2f/%.2f prop %.2f/%.2f adp %.2f/%.2f mpc %.2f/%.2f mpc* %.2f/%.2f\n', O, P, A, Q, Qs);
