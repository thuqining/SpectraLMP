% FIG_ROBUST3  Fig. 2, IEEE 14-bus system, extra annual cost over the
% hindsight optimum (% of it): (a) six fleets, 60-day history, bars without
% line limit and diamonds with line 5-6 at 3.1 kW; (b) history length,
% average over the six fleets and both line cases (MPC has 60 days only);
% (c) fleet 6w against the share of steps in which line 5-6 binds.
% Methods use red/grays so that they differ from the unit colors of Fig. 1.
% Writes Figure/fig_robust3.pdf.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
ld = @(t) getfield(load(fullfile(root, 'results', 'net', ['i14_' t '.mat'])), 'res'); %#ok<GFLD>
pc = @(r, f) 100*r.(f)/r.lb - 100;
ft = {'free', 'lim'};
% (a) fleets, 60-day history
V = nan(6, 3, 2);
for l = 1:2
    for f = 1:6
        r = ld(sprintf('fleet_%s_60_%d', ft{l}, f));
        if f == 4, q = ld(sprintf('mpcfit_%s', ft{l}));
        elseif l == 1, q = ld(sprintf('mpcfit_free_fleet_%d', f));
        else, q = ld(sprintf('mpcfit_fleet_%d', f)); end
        V(f,:,l) = [pc(r,'prop'), pc(r,'adp'), pc(q,'mpcfit')];
    end
end
% (b) history length: mean over fleets and line cases
hs = [7 30 60]; Hh = zeros(3, 2);
for k = 1:3
    v = zeros(12, 2); n = 0;
    for l = 1:2
        for f = 1:6
            r = ld(sprintf('fleet_%s_%d_%d', ft{l}, hs(k), f)); n = n + 1;
            v(n,:) = [pc(r,'prop'), pc(r,'adp')];
        end
    end
    Hh(k,:) = mean(v);
end
% (c) fleet 6w without limit and at 4.2, 3.1, 2.4, 1.5 kW
L = {ld('fleet_free_60_6'), ld('wide6_4.2_1'), ld('fleet_lim_60_6'), ld('wide6_2.4_1'), ld('wide6_1.5_1')};
Q = {ld('mpcfit_free_fleet_6'), ld('wide6_4.2_2'), ld('mpcfit_fleet_6'), ld('wide6_2.4_2'), ld('wide6_1.5_2')};
B = zeros(5, 3); sh = zeros(1, 5);
for k = 1:5
    B(k,:) = [pc(L{k},'prop'), pc(L{k},'adp'), pc(Q{k},'mpcfit')]; sh(k) = 100*L{k}.congested;
end

fn = 'Times New Roman'; fs = 8; yl = [0 40];
toRGB = @(h) sscanf(h(2:end), '%2x%2x%2x')' / 255;
col = cellfun(toRGB, {'#B2182B', '#4D4D4D', '#B3B3B3'}, 'UniformOutput', false);
W = 3.5; H = 1.10; in2n = @(l, b, w, h) [l/W, b/H, w/W, h/H];
f = figure('Units', 'inches', 'Position', [0.5 0.5 W H], 'Color', 'w', 'Visible', 'off');
opt = {'FontName', fn, 'FontSize', fs, 'LabelFontSizeMultiplier', 1, 'TickDir', 'out', 'Box', 'off', ...
    'TickLength', [0.02 0.02], 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 1};
% (a) y axis compressed above 40 (factor 1/3) so that the MPC bars of 4c
% and 6c and their congested diamonds keep their true order and position
g = @(y) min(y, yl(2)) + max(y - yl(2), 0)/3;
ax = axes(f, 'Position', in2n(0.33, 0.30, 1.62, 0.62)); hold(ax, 'on');
b = bar(ax, 1:6, g(V(:,:,1)), 'grouped', 'BarWidth', 0.84, 'EdgeColor', 'none');
for j = 1:3, b(j).FaceColor = col{j}; end
for j = 1:3
    xe = b(j).XEndPoints;
    hm = plot(ax, xe, g(V(:,j,2)), 'd', 'MarkerSize', 3.0, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'LineStyle', 'none');
    for i = 1:6
        if V(i,j,1) > yl(2)                              % break mark where the scale changes
            plot(ax, xe(i) + [-0.14 0.14], yl(2) - [1.0 -0.6], 'w-', 'LineWidth', 1.4);
        end
    end
end
set(ax, opt{:}, 'XLim', [0.5 6.5], 'XTick', 1:6, 'XTickLabel', {'2c', '2w', '4c', '4w', '6c', '6w'}, ...
    'YLim', [0 g(70)], 'YTick', [0:10:40 g(60)], 'YTickLabel', {'0', '10', '20', '30', '40', '60'});
% axis break mark on the y axis between 40 and 60
plot(ax, 0.5 + [-0.12 0.12], yl(2) + [0.4 1.6], 'k-', 'LineWidth', 0.6, 'Clipping', 'off');
plot(ax, 0.5 + [-0.12 0.12], yl(2) + [1.4 2.6], 'k-', 'LineWidth', 0.6, 'Clipping', 'off');
ylabel(ax, 'Extra cost (%)'); xlabel(ax, 'Fleet');
text(ax, 0.72, g(70) - 0.3, '(a)', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'top');
% (b)
ax2 = axes(f, 'Position', in2n(2.17, 0.30, 0.52, 0.62)); hold(ax2, 'on');
for j = 1:2
    plot(ax2, 1:3, Hh(:,j), '-o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 2.8, 'LineWidth', 1.0);
end
set(ax2, opt{:}, 'XLim', [0.6 3.4], 'XTick', 1:3, 'XTickLabel', {'7', '30', '60'}, 'YLim', [8 17], 'YTick', [8 12 16]);
xlabel(ax2, 'Initial days');
text(ax2, 0.7, 16.95, '(b)', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'top');
% (c)
ax3 = axes(f, 'Position', in2n(2.90, 0.30, 0.56, 0.62)); hold(ax3, 'on');
for j = 1:3
    plot(ax3, sh, B(:,j), '-o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 2.8, 'LineWidth', 1.0);
end
set(ax3, opt{:}, 'XLim', [-8 78], 'XTick', [0 35 70], 'YLim', [0 46], 'YTick', 0:20:40);
xlabel(ax3, 'Congestion (%)');
text(ax3, -6, 45.5, '(c)', 'FontName', fn, 'FontSize', fs, 'VerticalAlignment', 'top');
lg = legend(ax, [b, hm], {'Proposed', 'ADP', 'MPC', 'Congested'}, 'FontName', fn, 'FontSize', fs, ...
    'Box', 'off', 'Orientation', 'horizontal', 'NumColumns', 4);
lg.ItemTokenSize = [9 6]; lg.Position = in2n(0.33, 0.96, 3.13, 0.13);
exportgraphics(f, fullfile(root, 'Figure', 'fig_robust3.pdf'), 'ContentType', 'vector'); close(f);
disp(round(Hh, 2)); disp(round([sh' B], 1));
