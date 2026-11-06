% RUN_HINDCAL  Consistent hindsight oracle for Table I (reviewer R1-C1):
% the proposed method with its own past-fitted calibration and fallback,
% but with the predictions replaced by the exact test-year averages.
% Writes results/net/i14_hindcal_free.mat and i14_hindcal_lim.mat.
for k = 1:2
    cfg = netcfg(k, 'main14');
    cfg.methods = {'prop', 'hindcal'};
    cfg.tag = strrep(cfg.tag, 'main', 'hindcal');
    t0 = tic; res = netrun(cfg);
    fprintf('[run_hindcal] %s: prop %.2f%%, hindcal %.2f%% (%.0f s)\n', cfg.tag, ...
        100*res.prop/res.lb, 100*res.hindcal/res.lb, toc(t0));
end
