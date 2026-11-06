% RUN_MPC_HIST  Rerun every fitted-forecast MPC result of the letter with the
% forecast lags taken from the last training day instead of the test-period
% mean (reviewer R1, minor 8). Saves results/net/<tag>_h.mat next to the
% original files and prints old versus new costs (% of hindsight optimum).
root = fileparts(fileparts(mfilename('fullpath')));
ld = @(t) getfield(load(fullfile(root, 'results', 'net', [t '.mat'])), 'res'); %#ok<GFLD>
jobs = {};
for id = 1:2, jobs{end+1} = netcfg(id, 'mpcfit14'); end                 %#ok<SAGROW> Table I
for id = [1 2 3 5 6], jobs{end+1} = netcfg(id, 'mpcfitfig14'); end      %#ok<SAGROW> Fig. 2(a), limited
for id = [1 2 3 5 6], jobs{end+1} = netcfg(id, 'mpcfitfree14'); end     %#ok<SAGROW> Fig. 2(a), free
for id = 4:6, jobs{end+1} = netcfg(id, 'wide614'); end                  %#ok<SAGROW> Fig. 2(c)
for k = 1:numel(jobs)
    cfg = jobs{k}; old = cfg.tag; cfg.tag = [old '_h'];
    t0 = tic; r = netrun(cfg);
    o = ld(old);
    fprintf('[mpc_hist] %-24s old %.2f%%  new %.2f%%  (%.0f s)\n', old, 100*o.mpcfit/o.lb, 100*r.mpcfit/r.lb, toc(t0));
end
