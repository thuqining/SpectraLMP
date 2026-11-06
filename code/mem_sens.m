% MEM_SENS  Memory length of the predictors (hours of past net load) on the
% uncongested IEEE 14-bus case: cost of the proposed policy, % of hindsight.
for mh = [48 96 384 768]
    cfg = netcfg(1, 'main14'); cfg.methods = {'prop'}; cfg.memH = mh;
    cfg.tag = sprintf('i14_mem%d_free', mh);
    netrun(cfg);
end
