% CHECK_SPECTRAL_SHARE  Empirical frequency split of the proposed dispatch
% (base fleet, no congestion, test year) against the shares H_i of (5).
root = fileparts(fileparts(mfilename('fullpath')));
S = load(fullfile(root, 'results', 'net', 'i14_frame_free.mat')); tr = S.res.traj;
disp(fieldnames(tr)');
x = sum(tr.Xtest, 2); P = tr.p; if size(P,1) ~= numel(x), P = P'; end
fprintf('T = %d, N = %d\n', numel(x), size(P, 2));
G = load(fullfile(root, 'results', 'net14_gap.mat')); pol = G.pol;
c = pol.ceff(:); k = pol.kappa(:); q = 0.0073; dt = 0.25;
% Welch cross-spectra, 30-day segments
nw = 96*30; win = hann(nw); nov = nw/2;
[Sxx, f] = cpsd(x - mean(x), x - mean(x), win, nov, nw, 1/dt);
cyc = f*24;                                   % cycles/day
for i = 1:size(P, 2)
    Spx = cpsd(P(:,i) - mean(P(:,i)), x - mean(x), win, nov, nw, 1/dt);
    He(:,i) = real(Spx ./ Sxx);                %#ok<SAGROW>
end
w = 2*pi*f; a = c' + k'./max(w,1e-9).^2; Ht = (1./a) ./ (1/q + sum(1./a, 2));
sel = cyc >= 0.05 & cyc <= 40;
for b = [0.1 0.3 1 3 10 30]
    [~, j] = min(abs(cyc - b));
    fprintf('%5.1f c/d | emp %s | thy %s\n', cyc(j), mat2str(round(He(j,:),2)), mat2str(round(Ht(j,:),2)));
end
save(fullfile(root, 'results', 'spectral_share.mat'), 'cyc', 'He', 'Ht', 'sel');
