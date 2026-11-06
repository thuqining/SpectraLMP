function [X, Lv] = net14_loads(Xb, net, nTrain, Xraw)
% NET14_LOADS  Net-load deviations at the eleven load buses of net14_setup.
% Three shapes: the total of the three buildings (buses 2, 3, 4),
% industrial 1 (buses 5, 6, 9, 10), and residential 6 (buses 11-14). Each
% shape is scaled to the training deviation of the total, and bus n takes
% the share Pd_n/sum(Pd) of its shape, so bus magnitudes keep the case14
% load ratios.
%   Xb    T x 3 building deviations [ind1 res5 res6] (each minus its causal
%         weekly mean); nTrain training samples (scaling uses them only).
%   Xraw  (optional) T x 3 building net loads; Lv returns the bus net
%         loads with the same shapes, scales, and shares (for forecast
%         errors relative to the net load itself).
sh = [sum(Xb, 2), Xb(:,1), Xb(:,3)];
sx = std(sh(1:nTrain, 1));
sc = sx ./ std(sh(1:nTrain, :));
sh = sh .* sc;
grp = containers.Map({2,3,4, 5,6,9,10, 11,12,13,14}, {1,1,1, 2,2,2,2, 3,3,3,3});
w = net.Pd / sum(net.Pd);
X = zeros(size(Xb, 1), numel(net.loadBus)); Lv = [];
if nargin > 3, shr = [sum(Xraw, 2), Xraw(:,1), Xraw(:,3)] .* sc; Lv = X; end
for k = 1:numel(net.loadBus)
    X(:, k) = w(k) * sh(:, grp(net.loadBus(k)));
    if nargin > 3, Lv(:, k) = w(k) * shr(:, grp(net.loadBus(k))); end
end
end
