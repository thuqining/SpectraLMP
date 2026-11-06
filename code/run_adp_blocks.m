function out = run_adp_blocks(x, par, fits, blk, name)
% RUN_ADP_BLOCKS  ADP dispatch with values refitted per block (fits{b}).
if nargin < 5, name = 'ADP (separable values)'; end
T = numel(x); p = zeros(par.N,T); s = zeros(1,T); e = zeros(par.N,T); ep = par.e0(:)';
tic;
for t = 1:T
    f = fits{ceil(t/blk)};
    [pt, st, ep] = dispatch_value_step(x(t), ep, max(f(1,:) + f(2,:).*ep, 0), par);
    p(:,t) = pt'; s(t) = st; e(:,t) = ep';
end
tt = toc;
[~, out] = evalc('pack_results(name, p, s, e, x, par, tt/T*1e3)');
end
