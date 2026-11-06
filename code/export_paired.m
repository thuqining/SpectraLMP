% EXPORT_PAIRED  Paired result tables for the SpectraLMP release:
% paired_36.csv (six fleets x two line cases x three histories) and
% paired_24.csv (base case and 11 one-at-a-time changes x two line cases).
% Costs are in % of each case's hindsight optimum.
root = fileparts(fileparts(mfilename('fullpath')));
ld = @(t) getfield(load(fullfile(root, 'results', 'net', ['i14_' t '.mat'])), 'res'); %#ok<GFLD>
pc = @(r, f) 100*r.(f)/r.lb;
out = fullfile(root, 'SpectraLMP', 'results'); if ~isfolder(fullfile(root, 'SpectraLMP')), out = fullfile(root, 'results'); end
if ~isfolder(out), mkdir(out); end
names = {'2c', '2w', '4c', '4w', '6c', '6w'}; ft = {'free', 'lim'}; lines = {'none', '3.1 kW'};
fid = fopen(fullfile(out, 'paired_36.csv'), 'w');
fprintf(fid, 'fleet,line_limit,history_days,proposed,adp,mpc,adp_minus_proposed,mpc_minus_proposed\n');
for l = 1:2
    for h = [60 30 7]
        for f = 1:6
            r = ld(sprintf('fleet_%s_%d_%d', ft{l}, h, f));
            m = NaN;
            if h == 60
                if f == 4, q = ld(sprintf('mpcfit_%s', ft{l}));
                elseif l == 1, q = ld(sprintf('mpcfit_free_fleet_%d', f));
                else, q = ld(sprintf('mpcfit_fleet_%d', f)); end
                m = pc(q, 'mpcfit');
            end
            fprintf(fid, '%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n', names{f}, lines{l}, h, pc(r,'prop'), pc(r,'adp'), ...
                m, pc(r,'adp') - pc(r,'prop'), m - pc(r,'prop'));
        end
    end
end
fclose(fid);
fid = fopen(fullfile(out, 'paired_24.csv'), 'w');
fprintf(fid, 'setting,line_limit,proposed,adp,adp_minus_proposed\n');
for id = 1:24
    cfg = netcfg(id, 'sens14'); r = ld(strrep(cfg.tag, 'i14_', ''));
    lab = cfg.label; if isfinite(cfg.F), ln = '3.1 kW'; else, ln = 'none'; end
    fprintf(fid, '%s,%s,%.2f,%.2f,%.2f\n', lab, ln, pc(r,'prop'), pc(r,'adp'), pc(r,'adp') - pc(r,'prop'));
end
fclose(fid);
disp('exported');
