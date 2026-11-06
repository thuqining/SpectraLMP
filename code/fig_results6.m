% FIG_RESULTS6  Fig. 2, IEEE 14-bus system, three panels (3.5 x 2.96 in):
% (a) hindsight LMP without congestion (gray, $/MWh) and the average (4)
%     of each unit (color) over a window matched to its time constant
%     (2 days, 1 week, 2 months, 1 year, nested around the week of
%     Fig. 1); red where the LMP is above (discharge), blue below (charge);
% (b) annual cost in % of the hindsight optimum of the proposed method,
%     ADP, and MPC in twelve tests (six fleets with 60 days of history,
%     without and with the line limit); MPC bars beyond the axis are
%     capped and labeled;
% (c) storage share of each unit by frequency: theory H_i/sum H of (5)
%     (lines) against the dispatch of Algorithm 1 over the test year
%     (dots), with the crossovers.
% Writes Figure/fig_results6.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
fn = 'Times New Roman'; fs = 7.5;
toRGB = @(h) sscanf(h(2:end), '%2x%2x%2x')' / 255;
cP = toRGB('#B2182B'); cA = toRGB('#2166AC'); cM = toRGB('#B8B8B8');
cu = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;      % unit colors of Fig. 1
W = 3.5; H = 2.45; in2n = @(l, b, w, h) [l/W, b/H, w/W, h/H];
f = figure('Units', 'inches', 'Position', [0.5 0.5 W H], 'Color', 'w', 'Visible', 'off');
opt = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'TickDir', 'out', 'Box', 'off', ...
    'TickLength', [0.012 0.012]};
ul = {'1 h', '4 h', '30 h', '200 h'};

% ---- (a) LMP and the average of each unit ----
D = load(fullfile(root, 'results', 'fig_idea_s.mat'), 'lam', 'tau', 'dt'); spd = 96;
lam = 1000*D.lam(:); T = numel(lam);                % $/MWh
fq = (0:T-1)'/T; fq(fq > 0.5) = fq(fq > 0.5) - 1; wq = 2*pi*fq/D.dt;
mu = zeros(T, 4);
for i = 1:4, mu(:, i) = real(ifft(fft(lam) ./ (1 + D.tau(i)^2 * wq.^2))); end
win = 7*spd; best = 1; bv = -inf;                   % same week as Fig. 1
for s0 = 30*spd:spd:T - 30*spd - win
    v = std(mu(s0:s0+win-1, 3)); if v > bv, bv = v; best = s0; end
end
mid = best + win/2; len = [2 7 60 365]*spd; bin = [1 1 8 96];
xl = {{'0', '0.5', '1', '1.5', '2 d'}, {'0', '1', '2', '3', '4', '5', '6', '7 d'}, ...
      {'0', '10', '20', '30', '40', '50', '60 d'}, {'0', '2', '4', '6', '8', '10', '12 mo'}};
xk = cellfun(@(c) linspace(0, 1, numel(c)), xl, 'UniformOutput', false);
cD = [214 96 77]/255; cC = [67 147 195]/255;
ha = 0.14; ga = 0.142; yt0 = 2.27;
for i = 1:4
    s0 = max(1, min(T - len(i) + 1, round(mid - len(i)/2))); id = s0:s0+len(i)-1;
    nb = floor(len(i)/bin(i)); id = id(1:nb*bin(i));
    L = mean(reshape(lam(id), bin(i), nb), 1)'; M = mean(reshape(mu(id, i), bin(i), nb), 1)';
    tt = ((1:nb)' - 0.5)/nb;
    axa = axes(f, 'Position', in2n(0.42, yt0 - i*ha - (i - 1)*ga, 2.78, ha)); hold(axa, 'on');
    up = max(L, M); dn = min(L, M);
    fill(axa, [tt; flipud(tt)], [up; flipud(M)], cD, 'EdgeColor', 'none', 'FaceAlpha', 0.40);
    fill(axa, [tt; flipud(tt)], [M; flipud(dn)], cC, 'EdgeColor', 'none', 'FaceAlpha', 0.40);
    plot(axa, tt, M, '-', 'Color', cu(i,:), 'LineWidth', 1.1);          % average under the LMP
    plot(axa, tt, L, '-', 'Color', [0 0 0], 'LineWidth', 0.45);
    yr = [min(L) max(L)]; ytk = ceil(yr(1)/100)*100:100:floor(yr(2)/100)*100;
    if numel(ytk) > 2, ytk = ytk([1 end]); end
    set(axa, 'FontName', fn, 'FontSize', fs - 2, 'XLim', [0 1], 'YLim', yr + [-0.06 0.06]*diff(yr), ...
        'XTick', xk{i}, 'XTickLabel', xl{i}, 'YTick', ytk, 'TickDir', 'out', 'TickLength', [0.006 0.006], ...
        'Box', 'on', 'XColor', [0.35 0.35 0.35], 'YColor', [0.35 0.35 0.35], 'Color', 'w');
    text(axa, 1.015, 0.5, ul{i}, 'Units', 'normalized', 'Color', cu(i,:), 'FontName', fn, ...
        'FontSize', fs - 0.5, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
xlabel(axa, 'Time', 'FontName', fn, 'FontSize', fs, 'Color', [0 0 0]);
axa.XLabel.Units = 'normalized'; axa.XLabel.Position(2) = -0.62;

% ---- (b) cost by test case ----
rel = fullfile(root, 'results');                      % paired_36.csv (export_paired_s)
Tb = readtable(fullfile(rel, 'paired_36.csv'), 'TextType', 'string');
Tb = Tb(Tb.history_days == 60, :);
fl = ["2c" "2w" "4c" "4w" "6c" "6w"]; lns = ["none" "3.3 kW"]; cap = 127; bw = 0.25;
ax3 = axes(f, 'Position', in2n(0.36, 0.34, 1.40, 0.56)); hold(ax3, 'on');
bx = @(x, y1, c) patch(ax3, x + bw/2*[-1 1 1 -1], [100 100 y1 y1], c, 'EdgeColor', 'none');
mk = @(x, y1) plot(ax3, x, y1, 'o', 'MarkerSize', 2.6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.6);
xt = 1:6; vv = @(r, m) r.(m);
for k = 1:6
    r1 = Tb(Tb.fleet == fl(k) & Tb.line_limit == lns(1), :);    % bars: no line limit
    r2 = Tb(Tb.fleet == fl(k) & Tb.line_limit == lns(2), :);    % marks: congested
    m3 = {'proposed', 'adp', 'mpc'}; cc = {cP, cA, cM};
    for j = 1:3
        x = k + (j - 2)*bw; y1 = vv(r1, m3{j}); y2 = vv(r2, m3{j});
        bx(x, min(y1, cap), cc{j}); mk(x, min(y2, cap));
        if y1 > cap
            text(ax3, x, cap + 0.4, sprintf('%.0f', y1), 'FontName', fn, 'FontSize', fs - 2, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', [0.35 0.35 0.35]);
        end
    end
end
set(ax3, opt{:}, 'XLim', [0.5 6.5], 'XTick', xt, 'XTickLabel', cellstr(fl), 'YLim', [100 131], ...
    'YTick', 100:10:120, 'YGrid', 'on', 'GridColor', [0.9 0.9 0.9], 'GridAlpha', 1, 'Layer', 'top');
ax3.XAxis.FontSize = fs - 1; ax3.XAxis.TickLength = [0 0];
ylabel(ax3, 'Cost (%)'); xlabel(ax3, 'Fleet (c: 4–30 h, w: 1–200 h)');
ax3.XLabel.Units = 'normalized'; ax3.XLabel.Position(2) = -0.21; ax3.XLabel.FontSize = fs;
hP = patch(ax3, nan, nan, cP, 'EdgeColor', 'none'); hA = patch(ax3, nan, nan, cA, 'EdgeColor', 'none');
hM = patch(ax3, nan, nan, cM, 'EdgeColor', 'none'); hK = plot(ax3, nan, nan, 'o', 'MarkerSize', 2.6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.6);
lg = legend(ax3, [hP hA hM hK], {'Proposed', 'ADP', 'MPC', 'Congested'}, 'FontName', fn, ...
    'FontSize', fs - 1, 'Box', 'off', 'Orientation', 'horizontal', 'NumColumns', 2);
lg.ItemTokenSize = [7 6]; lg.Position = in2n(0.42, 0.91, 1.30, 0.22);

% ---- (c) storage share by frequency ----
S = load(fullfile(root, 'results', 'spectral_share_s.mat'));
G = load(fullfile(root, 'results', 'pol_s.mat')); c = G.pol.ceff(:)'; kap = G.pol.kappa(:)';
cy = logspace(log10(0.02), log10(48), 400)'; wv = 2*pi*cy/24;     % rad/h
Ha = (1 ./ (c + kap ./ wv.^2)); Hn = Ha ./ sum(Ha, 2);
He = max(S.He, 0); He = He ./ sum(He, 2);
ax2 = axes(f, 'Position', in2n(2.17, 0.34, 1.28, 0.56)); hold(ax2, 'on');
[~, dom] = max(Hn, [], 2); xc = [];
for k = 2:numel(cy), if dom(k) ~= dom(k-1), xc(end+1) = sqrt(cy(k)*cy(k-1)); end, end %#ok<AGROW>
for k = 1:numel(xc)
    plot(ax2, xc(k)*[1 1], [0 1], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    P = 24/xc(k); if P > 48, ps = sprintf('%.0f d', P/24); else, ps = sprintf('%.2g h', P); end
    text(ax2, xc(k)*1.08, 0.99, ps, 'FontName', fn, 'FontSize', fs - 2, 'Color', [0.35 0.35 0.35], ...
        'VerticalAlignment', 'top');
end
for i = 1:4, semilogx(ax2, cy, Hn(:, i), '-', 'Color', cu(i,:), 'LineWidth', 1.0); end
edges = logspace(log10(0.03), log10(40), 15); cb = sqrt(edges(1:end-1).*edges(2:end));
for i = 1:4
    hb = nan(size(cb));
    for k = 1:numel(cb)
        m = S.cyc >= edges(k) & S.cyc < edges(k+1); if any(m), hb(k) = mean(He(m, i)); end
    end
    plot(ax2, cb, hb, 'o', 'MarkerSize', 2.2, 'MarkerFaceColor', cu(i,:), 'MarkerEdgeColor', 'none');
end
set(ax2, opt{:}, 'XScale', 'log', 'XLim', [0.02 48], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'}, ...
    'YLim', [0 1], 'YTick', [0 0.5 1], 'Layer', 'top');
ylabel(ax2, 'Share'); xlabel(ax2, 'Frequency (cycles/day)');
ax2.XLabel.Units = 'normalized'; ax2.XLabel.Position(2) = -0.21;
lc = [12 1.8 0.13 0.028]; dy = [0.05 0.05 0.05 -0.40];
for i = 1:4
    [~, k] = min(abs(cy - lc(i)));
    text(ax2, lc(i), Hn(k, i) + dy(i), ul{i}, 'Color', cu(i,:), 'FontName', fn, 'FontSize', fs - 1, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
hT = plot(ax2, nan, nan, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
hE = plot(ax2, nan, nan, 'o', 'MarkerSize', 2.6, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none');
lg2 = legend(ax2, [hT hE], {'Hindsight', 'Proposed'}, 'FontName', fn, 'FontSize', fs - 1, 'Box', 'off', ...
    'Orientation', 'horizontal', 'NumColumns', 2);
lg2.ItemTokenSize = [10 6]; lg2.Position = in2n(2.17, 0.95, 1.28, 0.11);

% ---- (a) key, panel letters, and the LMP axis label ----
axo = axes(f, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 W], 'YLim', [0 H]); hold(axo, 'on');
yk = 2.36;
plot(axo, [0.55 0.75], [yk yk], '-', 'Color', [0 0 0], 'LineWidth', 0.45);
text(axo, 0.78, yk, 'LMP', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
plot(axo, [1.12 1.32], [yk yk], '-', 'Color', cu(1,:), 'LineWidth', 1.1);
text(axo, 1.35, yk, 'LMP-average', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
patch(axo, [2.05 2.17 2.17 2.05], yk + [-0.04 -0.04 0.04 0.04], cD, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
text(axo, 2.20, yk, 'Discharge', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
patch(axo, [2.75 2.87 2.87 2.75], yk + [-0.04 -0.04 0.04 0.04], cC, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
text(axo, 2.90, yk, 'Charge', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
text(axo, 0.07, (yt0 + yt0 - 4*ha - 3*ga)/2, 'LMP ($/MWh)', 'FontName', fn, 'FontSize', fs, 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
pl = {'(a)', 0.02, yk; '(b)', 0.02, 1.06; '(c)', 1.86, 1.06};
for k = 1:3
    text(axo, pl{k,2}, pl{k,3}, pl{k,1}, 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
end
exportgraphics(f, fullfile(root, 'Figure', 'fig_results6.pdf'), 'ContentType', 'vector'); close(f);
fprintf('crossover periods (h): %s\n', mat2str(24./xc, 3));
