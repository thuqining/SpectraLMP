% PICK_WEEK  Score each test week by how distinct the four SoC curves are:
% the fast units should cycle, the 30-h unit should move slowly, and the
% 200-h unit should hold energy. Prints the best candidates.
root = fileparts(fileparts(mfilename('fullpath')));
R = load(fullfile(root, 'results', 'net', 'i14s_main_free.mat')); tr = R.res.traj;
soc = tr.e' ./ tr.E(:)'; T = size(soc, 1); spd = 96; win = 7*spd;
best = [];
for s = 1:spd:T - win
    w = soc(s:s+win-1, :);
    sd = std(w); mn = mean(w);
    daily = mean(abs(diff(w(1:spd:end, :))));          % day-to-day change
    sc = min(sd(1), 0.4) + min(sd(2), 0.4) + 2*min(range(w(:,3)), 0.5) + 2*min(mn(4), 0.4) - abs(sd(1) - sd(2));
    best(end+1, :) = [s, sc, sd, mn(4), range(w(:,3))]; %#ok<AGROW>
end
[~, o] = sort(best(:, 2), 'descend');
fprintf('day  score  sd1  sd2  sd3  sd4  mean4  range3\n');
for k = 1:8, b = best(o(k), :); fprintf('%5.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', (b(1)-1)/spd, b(2:end)); end
