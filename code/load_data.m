function [x, dt, meta] = load_data(csvPath, monthsNeeded, keepB)
% LOAD_DATA  Build a community-microgrid net-load series from OPSD CoSSMic data.
%   Cumulative meter readings (kWh) are differenced into 15-min average power
%   (kW). Rows whose value was gap-filled by OPSD (flagged in the
%   'interpolated' column) are treated as invalid for the columns concerned.
%   A greedy column-drop keeps the largest building set that still yields a
%   contiguous window of monthsNeeded. Community load = sum(grid_import + pv
%   - export); community PV = sum(pv), scaled so mean(net load) = 0.

dt = 0.25;
stepsPerDay = 96;
need = round(monthsNeeded * 30 * stepsPerDay);

% building -> columns map (buildings are dropped as a whole)
B = struct();
B.ind1 = {{'DE_KN_industrial1_grid_import'}, {'DE_KN_industrial1_pv_1','DE_KN_industrial1_pv_2'}, {}};
B.ind3 = {{'DE_KN_industrial3_grid_import'}, {'DE_KN_industrial3_pv_facade','DE_KN_industrial3_pv_roof'}, {}};
B.pub1 = {{'DE_KN_public1_grid_import'}, {}, {}};
B.pub2 = {{'DE_KN_public2_grid_import'}, {}, {}};
B.res1 = {{'DE_KN_residential1_grid_import'}, {'DE_KN_residential1_pv'}, {}};
B.res2 = {{'DE_KN_residential2_grid_import'}, {}, {}};
B.res3 = {{'DE_KN_residential3_grid_import'}, {'DE_KN_residential3_pv'}, {'DE_KN_residential3_grid_export'}};
B.res4 = {{'DE_KN_residential4_grid_import'}, {'DE_KN_residential4_pv'}, {'DE_KN_residential4_grid_export'}};
B.res5 = {{'DE_KN_residential5_grid_import'}, {}, {}};
B.res6 = {{'DE_KN_residential6_grid_import'}, {'DE_KN_residential6_pv'}, {'DE_KN_residential6_grid_export'}};
bn = fieldnames(B);

allCols = {};
for k = 1:numel(bn)
    for j = 1:3, allCols = [allCols, B.(bn{k}){j}]; end %#ok<AGROW>
end

opts = detectImportOptions(csvPath);
keep = [{'utc_timestamp','interpolated'}, allCols];
opts.SelectedVariableNames = intersect(opts.VariableNames, keep, 'stable');
opts = setvartype(opts, intersect(allCols, opts.VariableNames), 'double');
opts = setvartype(opts, {'interpolated'}, 'char');
T = readtable(csvPath, opts);
ts = T.utc_timestamp;
if ~isdatetime(ts), ts = datetime(ts, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ssX', 'TimeZone', 'UTC'); end
n = height(T) - 1;
interp = T.interpolated(2:end);          % flag applies to the step's end row

% per-column power and validity
P = struct(); V = struct();
for k = 1:numel(allCols)
    c = allCols{k};
    p = diff(T.(c)) / dt;
    bad = isnan(p) | p < 0;
    flagged = contains(interp, c);       % OPSD gap-filled -> not real dynamics
    p(bad) = NaN;
    p(flagged) = NaN;
    P.(c) = p;
    V.(c) = ~bad & ~flagged;
end

% validity per building
BV = false(n, numel(bn));
for k = 1:numel(bn)
    cols = [B.(bn{k}){:}];
    v = true(n,1);
    for j = 1:numel(cols), v = v & V.(cols{j}); end
    BV(:,k) = v;
end

% fixed building set (keepB), or greedy drop: keep as many buildings as
% possible subject to window >= need
active = true(1, numel(bn));
if nargin >= 3 && ~isempty(keepB)
    active = ismember(bn', keepB);
    assert(sum(active) == numel(keepB), 'Unknown building in keepB.');
    len = longest_run(all(BV(:,active), 2));
    assert(len >= need, 'Window of %d steps < %d needed.', len, need);
end
minPV = 2;                                % keep at least 2 PV buildings
hasPV = cellfun(@(f) ~isempty(B.(f){2}), bn)';
while true
    [len, i0] = longest_run(all(BV(:,active), 2));
    if len >= need, break; end
    % try dropping each active building, take the drop with max window
    best = 0; bestk = 0;
    for k = find(active)
        if hasPV(k) && sum(hasPV & active) <= minPV, continue; end
        if sum(active) <= 4, continue; end
        tryA = active; tryA(k) = false;
        l2 = longest_run(all(BV(:,tryA), 2));
        if l2 > best, best = l2; bestk = k; end
    end
    assert(bestk > 0, 'Cannot reach %d contiguous steps (best %d).', need, len);
    active(bestk) = false;
    fprintf('[data] dropping %s (window -> %d steps)\n', bn{bestk}, best);
end
[~, i0] = longest_run(all(BV(:,active), 2));
idx = i0 : i0 + need - 1;

L = zeros(need,1); G = zeros(need,1);
act = find(active); Lb = zeros(need, numel(act)); Gb = Lb;
for kk = 1:numel(act)
    k = act(kk);
    lc = B.(bn{k}){1}; pc = B.(bn{k}){2}; ec = B.(bn{k}){3};
    for j = 1:numel(lc), Lb(:,kk) = Lb(:,kk) + P.(lc{j})(idx); end
    for j = 1:numel(pc), Lb(:,kk) = Lb(:,kk) + P.(pc{j})(idx); Gb(:,kk) = Gb(:,kk) + P.(pc{j})(idx); end
    for j = 1:numel(ec), Lb(:,kk) = Lb(:,kk) - P.(ec{j})(idx); end
    L = L + Lb(:,kk); G = G + Gb(:,kk);
end
alpha = mean(L) / mean(G);
x = L - alpha * G;
xb = Lb - alpha * Gb;                     % per-building net load, sum(xb,2) = x
assert(std(x) > 1, 'Degenerate net load (std = %.3f kW).', std(x));

meta = struct('t0', ts(idx(1)+1), 't1', ts(idx(end)+1), ...
    'meanLoad', mean(L), 'stdX', std(x), 'alpha', alpha, ...
    'nSteps', numel(x), 'buildings', {bn(active)}, 'xb', xb);
fprintf('[data] %s -> %s | %d buildings | mean load %.1f kW | PV scale %.2f | std(x) %.1f kW\n', ...
    string(meta.t0), string(meta.t1), sum(active), meta.meanLoad, alpha, meta.stdX);
end

function [len, i0] = longest_run(v)
dv = diff([0; v(:); 0]);
rs = find(dv == 1); re = find(dv == -1) - 1;
if isempty(rs), len = 0; i0 = 1; return; end
[len, k] = max(re - rs + 1);
i0 = rs(k);
end
