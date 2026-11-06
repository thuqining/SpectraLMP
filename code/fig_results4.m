% FIG_RESULTS4  Fig. 2, IEEE 14-bus system, three panels (3.5 x 2.75 in):
% (a) hindsight LMP without congestion (gray) and the average (4) of each
%     unit (color) over a window matched to its time constant (2 days,
%     1 week, 2 months, 1 year, nested around the week of Fig. 1); red
%     where the LMP is above (discharge), blue below (charge);
% (b) extra cost over the hindsight optimum of six fleets with 60 days of
%     history: proposed, ADP, and MPC (values beyond the axis in text);
%     for fleet 4w also the oracle with exact averages and MPC with
%     perfect forecasts; circles without congestion, diamonds with it;
% (c) storage share of each unit by frequency: theory H_i/sum H (lines)
%     against the dispatch of Algorithm 1 over the test year (dots), with
%     the band where each unit has the largest share and the crossovers.
% Writes Figure/fig_results4.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
ld = @(t) getfield(load(fullfile(root, 'results', 'net', ['i14_' t '.mat'])), 'res'); %#ok<GFLD>
pc = @(r, f) 100*r.(f)/r.lb - 100;
fn = 'Times New Roman'; fs = 7.5;
toRGB = @(h) sscanf(h(2:end), '%2x%2x%2x')' / 255;
cP = toRGB('#B2182B'); cA = toRGB('#4D4D4D'); cM = toRGB('#9E9E9E');
cu = [0 114 178; 230 159 0; 0 158 115; 204 121 167] / 255;      % unit colors of Fig. 1
W = 3.5; H = 2.75; in2n = @(l, b, w, h) [l/W, b/H, w/W, h/H];
f = figure('Units', 'inches', 'Position', [0.5 0.5 W H], 'Color', 'w', 'Visible', 'off');
opt = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'TickDir', 'out', 'Box', 'off', ...
    'TickLength', [0.012 0.012]};

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
ha = 0.145; ga = 0.02; ya0 = 1.98;
for i = 1:4
    s0 = max(1, min(T - len(i) + 1, round(mid - len(i)/2))); id = s0:s0+len(i)-1;
    nb = floor(len(i)/bin(i)); id = id(1:nb*bin(i));
    L = mean(reshape(lam(id), bin(i), nb), 1)'; M = mean(reshape(mu(id, i), bin(i), nb), 1)';
    tt = ((1:nb)' - 0.5)/nb;
    axa = axes(f, 'Position', in2n(0.60, ya0 + (4 - i)*(ha + ga), 2.85, ha)); hold(axa, 'on');
    up = max(L, M); dn = min(L, M);
    fill(axa, [tt; flipud(tt)], [up; flipud(M)], cD, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    fill(axa, [tt; flipud(tt)], [M; flipud(dn)], cC, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    plot(axa, tt, L, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.4);
    plot(axa, tt, M, '-', 'Color', cu(i,:), 'LineWidth', 0.9);
    yr = [min(L) max(L)];
    set(axa, 'XLim', [0 1], 'YLim', yr + [-0.06 0.06]*diff(yr), 'XTick', [], 'YTick', [], 'Box', 'on', ...
        'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], 'Color', 'w');
    text(axa, -0.115, 0.5, ul{i}, 'Units', 'normalized', 'Color', cu(i,:), 'FontName', fn, ...
        'FontSize', fs - 0.5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    text(axa, -0.012, 0.5, wl{i}, 'Units', 'normalized', 'Color', [0.35 0.35 0.35], 'FontName', fn, ...
        'FontSize', fs - 0.5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

% ---- (b) extra cost of six fleets ----
Mf = ld('main_free'); Ml = ld('main_lim'); Of = ld('hindcal_free'); Ol = ld('hindcal_lim');
O = [pc(Of,'hindcal'), pc(Ol,'hindcal')]; Qs = [pc(Mf,'mpcexact'), pc(Ml,'mpcexact')];
rel = fullfile(root, 'SpectraLMP', 'results'); if ~isfolder(rel), rel = fullfile(root, 'results'); end
Tb = readtable(fullfile(rel, 'paired_36.csv'), 'TextType', 'string');
Tb = Tb(Tb.history_days == 60, :);
fl = ["2c" "2w" "4c" "4w" "6c" "6w"]; lns = ["none" "3.1 kW"]; cap = 21; ms = 3.0;
ax3 = axes(f, 'Position', in2n(0.36, 1.12, 2.62, 0.68)); hold(ax3, 'on');
for k = 1:6
    mm = nan(1, 2);
    for l = 1:2
        r = Tb(Tb.fleet == fl(k) & Tb.line_limit == lns(l), :);
        y = k + (l - 1.5)*0.36; p = r.proposed - 100; a = r.adp - 100; mm(l) = r.mpc - 100;
        mk = 'o'; if l == 2, mk = 'd'; end
        if fl(k) == "4w"
            plot(ax3, [O(l) p], [y y], ':', 'Color', cP, 'LineWidth', 0.8);
            plot(ax3, O(l), y, mk, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cP, 'LineWidth', 0.7);
            plot(ax3, Qs(l), y, mk, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cM, 'LineWidth', 0.7);
        end
        plot(ax3, [p a], [y y], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.1);
        plot(ax3, a, y, mk, 'MarkerSize', ms, 'MarkerFaceColor', cA, 'MarkerEdgeColor', cA);
        plot(ax3, p, y, mk, 'MarkerSize', ms, 'MarkerFaceColor', cP, 'MarkerEdgeColor', cP);
        if mm(l) <= cap
            plot(ax3, mm(l), y, mk, 'MarkerSize', ms, 'MarkerFaceColor', cM, 'MarkerEdgeColor', cM);
        end
    end
    text(ax3, cap + 0.3, k, sprintf('%.0f/%.0f', mm), 'FontName', fn, 'FontSize', fs - 1, ...
        'VerticalAlignment', 'middle', 'Color', [0.35 0.35 0.35]);
end
text(ax3, cap + 0.3, 0.05, 'MPC', 'FontName', fn, 'FontSize', fs - 1, 'VerticalAlignment', 'middle', ...
    'Color', [0.35 0.35 0.35]);
set(ax3, opt{:}, 'YLim', [0.45 6.55], 'YTick', 1:6, 'YTickLabel', cellstr(fl), 'YDir', 'reverse', ...
    'XLim', [0 cap], 'XTick', 0:5:20, 'XGrid', 'on', 'GridColor', [0.9 0.9 0.9], 'GridAlpha', 1, ...
    'Clipping', 'off');
ax3.YAxis.FontSize = fs - 1;
xlabel(ax3, 'Extra cost (%)');
ax3.XLabel.Units = 'normalized'; ax3.XLabel.Position(2) = -0.20;
% legend
hP = plot(ax3, nan, nan, 'o', 'MarkerFaceColor', cP, 'MarkerEdgeColor', cP, 'MarkerSize', ms);
hA = plot(ax3, nan, nan, 'o', 'MarkerFaceColor', cA, 'MarkerEdgeColor', cA, 'MarkerSize', ms);
hM = plot(ax3, nan, nan, 'o', 'MarkerFaceColor', cM, 'MarkerEdgeColor', cM, 'MarkerSize', ms);
hO = plot(ax3, nan, nan, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cP, 'MarkerSize', ms, 'LineWidth', 0.7);
hS = plot(ax3, nan, nan, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cM, 'MarkerSize', ms, 'LineWidth', 0.7);
lg = legend(ax3, [hP hA hM hO hS], {'Proposed', 'ADP', 'MPC', 'Oracle', 'MPC*'}, 'FontName', fn, ...
    'FontSize', fs, 'Box', 'off', 'Orientation', 'horizontal', 'NumColumns', 5);
lg.ItemTokenSize = [7 6]; lg.Position = in2n(0.42, 2.62, 3.04, 0.13);

% ---- (c) storage share by frequency ----
S = load(fullfile(root, 'results', 'spectral_share.mat'));
G = load(fullfile(root, 'results', 'net14_gap.mat')); c = G.pol.ceff(:)'; kap = G.pol.kappa(:)'; q = 0.0073;
cy = logspace(log10(0.02), log10(48), 400)'; wv = 2*pi*cy/24;     % rad/h
Ha = (1 ./ (c + kap ./ wv.^2)); Hn = Ha ./ sum(Ha, 2);
He = max(S.He, 0); He = He ./ sum(He, 2);
ax2 = axes(f, 'Position', in2n(0.40, 0.30, 3.00, 0.52)); hold(ax2, 'on');
[~, dom] = max(Hn, [], 2); xc = [];
for k = 2:numel(cy), if dom(k) ~= dom(k-1), xc(end+1) = sqrt(cy(k)*cy(k-1)); end, end %#ok<AGROW>
ed = [cy(1) xc cy(end)]; dl = dom([1 find(diff(dom))' + 1]);
for k = 1:numel(dl)
    patch(ax2, [ed(k) ed(k+1) ed(k+1) ed(k)], [0 0 1 1], cu(dl(k),:), 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
for k = 1:numel(xc)
    plot(ax2, xc(k)*[1 1], [0 1], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    P = 24/xc(k); if P > 48, ps = sprintf('%.0f d', P/24); else, ps = sprintf('%.2g h', P); end
    text(ax2, xc(k)*1.08, 0.97, ps, 'FontName', fn, 'FontSize', fs - 1.5, 'Color', [0.35 0.35 0.35], ...
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
ax2.XLabel.Units = 'normalized'; ax2.XLabel.Position(2) = -0.30;
% panel letters
axo = axes(f, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 W], 'YLim', [0 H]); hold(axo, 'on');
pl = {'(a)', 0.02, 2.685; '(b)', 0.02, 1.90; '(c)', 0.02, 0.90};
for k = 1:3
    text(axo, pl{k,2}, pl{k,3}, pl{k,1}, 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'middle');
end
exportgraphics(f, fullfile(root, 'Figure', 'fig_results4.pdf'), 'ContentType', 'vector'); close(f);
fprintf('crossover periods (h): %s\n', mat2str(24./xc, 3));
