% SEAS_TABLE  Table I numbers of the seasonal 14-bus study (i14s_main_*).
R = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'net');
nm = {'free', 'lim'};
m = {'prop', 'adp', 'cut', 'mpcfit', 'rl', 'hindcal', 'mpcexact', 'mpc', 'myop', 'const', 'fixed', 'blind', 'common', 'hindavg'};
for k = 1:2
    S = load(fullfile(R, ['i14s_main_' nm{k} '.mat'])); r = S.res;
    fprintf('%s: congested %.1f%% | ', nm{k}, 100*r.congested);
    for j = 1:numel(m)
        if isfield(r, m{j}), fprintf('%s %.1f  ', m{j}, 100*r.(m{j})/r.lb); end
    end
    fprintf('\n   ms: ');
    f = fieldnames(r);
    for j = 1:numel(f), if startsWith(f{j}, 'ms'), fprintf('%s %.2f  ', f{j}, r.(f{j})); end, end
    if isfield(r, 'flowP'), fprintf('\n   flowP %.2f %.2f %.2f', r.flowP); end
    fprintf('\n');
end
