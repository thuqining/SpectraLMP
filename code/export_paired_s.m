% EXPORT_PAIRED_S  Seasonal 14-bus results (netcfg sets ending in '14s') for
% Fig. 2(b) and the text: paired_36.csv, six fleets x two line cases x three
% histories. proposed = Algorithm 1 (split predictor, set 'split14s');
% adp = ADP [4] with SoC-only values, adp_features = ADP with the features of
% Step 3, direct = hindsight duals as targets without the causal LMP average,
% averaged_target = the earlier design with the averages (4) as targets
% (set 'abl14s'); mpc = 24-h MPC with the XGBoost forecast (60 days only).
% Costs are in % of each case's hindsight optimum.
root = fileparts(fileparts(mfilename('fullpath')));
ld = @(t) getfield(load(fullfile(root, 'results', 'net', ['i14s_' t '.mat'])), 'res'); %#ok<GFLD>
pc = @(r, f) 100*r.(f)/r.lb;
out = fullfile(root, 'results');
if ~isfolder(out), mkdir(out); end
names = {'2c', '2w', '4c', '4w', '6c', '6w'}; ft = {'free', 'lim'}; lines = {'none', '3.3 kW'};
fid = fopen(fullfile(out, 'paired_36.csv'), 'w');
fprintf(fid, ['fleet,line_limit,history_days,proposed,adp,mpc,adp_minus_proposed,mpc_minus_proposed,' ...
    'adp_features,direct,averaged_target\n']);
for l = 1:2
    for h = [60 30 7]
        for f = 1:6
            s = ld(sprintf('split_%s_%d_%d', ft{l}, h, f));
            b = ld(sprintf('abl_%s_%d_%d', ft{l}, h, f));
            m = NaN;
            if h == 60
                q = ld(sprintf('mpcxgb_%s_fleet_%d', ft{l}, f));   % MPC with the XGBoost forecast
                m = pc(q, 'mpcxgb');
            end
            p = pc(s, 'split'); a = pc(b, 'adp');
            fprintf(fid, '%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', names{f}, lines{l}, h, p, a, ...
                m, a - p, m - p, pc(b, 'adpf'), pc(s, 'direct'), pc(b, 'prop'));
        end
    end
end
fclose(fid);
disp('exported');
