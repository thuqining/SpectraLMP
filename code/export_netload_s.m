% EXPORT_NETLOAD_S  Bus net loads of the seasonal 14-bus study (as netrun
% with trend = 0), saved to results/netload_s.mat for forecast_xgb.py.
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
Xs = meta.xb - mean(meta.xb, 1); nTrain = 60*spd;
net = net14_setup();
[X, ~] = net14_loads(Xs, net, nTrain, meta.xb);
hrs = mod(hour(meta.t0) + minute(meta.t0)/60 + (0:size(X,1)-1)'*dt, 24);
loadBus = net.loadBus;
save(fullfile(root, 'results', 'netload_s.mat'), 'X', 'nTrain', 'loadBus', 'hrs', 'dt');
fprintf('X %s, nTrain %d\n', mat2str(size(X)), nTrain);
