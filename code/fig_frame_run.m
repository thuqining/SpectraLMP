% FIG_FRAME_RUN  Proposed policy on the uncongested IEEE 14-bus case with
% trajectories saved (results/net/i14_frame_free.mat) for Fig. 1(a).
cfg = netcfg(1, 'main14'); cfg.methods = {'prop'}; cfg.tag = 'i14_frame_free'; cfg.traj = true;
netrun(cfg);
