function out = method_values_net(X, par, net, valueFun, name, opts)
% METHOD_VALUES_NET  Online network dispatch with stored-energy values.
% At each step, the operator solves a one-step DC-OPF QP in which unit i
% values its stored energy at v_i = valueFun(t, e) (step 2 of Algorithm 1
% under congestion): unit i discharges when its LMP exceeds
% v_i/etad_i + 2 c_i p_i and charges when it is below v_i etac_i - ...
%   X         T x numel(net.loadBus) nodal net loads (test period)
%   valueFun  @(t, e) -> 1 x N nonnegative values ($/kWh)
if nargin < 5, name = 'Online (network)'; end
% opts.emaTau (1 x N, h), opts.emaInit (1 x N): keep the causal exponential
% average Cm of the realized LMP at each unit's bus (up to step t-1) and call
% valueFun(t, e, Cm).
useEma = nargin >= 6 && isfield(opts, 'emaTau');
[T, ~] = size(X); N = par.N; dt = par.dt;
G = numel(net.genBus); C = numel(net.curtBus);
nv = 2*N + G + C;
lim = find(isfinite(net.F)); nL = numel(lim);
Gu = net.PTDF(lim, net.unitBus); Gg = net.PTDF(lim, net.genBus);
Gc = net.PTDF(lim, net.curtBus); Gx = net.PTDF(lim, net.loadBus);
c0 = par.c0(:)'; etac = par.etac(:)'; etad = par.etad(:)';
Q = sparse(blkdiag([diag(c0), diag(c0); diag(c0), diag(c0)], diag(net.genQ), net.qc*eye(C))) * dt;
Arow = [ones(1,N), -ones(1,N), ones(1,G), ones(1,C)];
Al = [Gu, -Gu, Gg, Gc];
model.Q = Q;
model.A = sparse([Arow; Al; Al]);
model.sense = [ '='; repmat('<', nL, 1); repmat('>', nL, 1) ];
params.OutputFlag = 0; params.Threads = 1;
vt = zeros(N,T); pd = zeros(N,T); pc = zeros(N,T); sg = zeros(G,T); sc = zeros(C,T); e = zeros(N,T);
ep = par.e0(:)'; lmp = zeros(T, net.nb);
if useEma, rE = exp(-dt ./ opts.emaTau(:)'); Cm = opts.emaInit(:)'; end
tic;
for t = 1:T
    if useEma, v = max(valueFun(t, ep, Cm), 0); else, v = max(valueFun(t, ep), 0); end
    vt(:,t) = v';
    pdmax = min(par.P, etad .* ep / dt);
    pcmax = min(par.P, (par.E - ep) ./ (etac * dt));
    model.obj = [ (v ./ etad) * dt, -(v .* etac) * dt, net.genA(:)' * dt, zeros(1, C) ]';
    xt = X(t,:);
    fx = Gx * xt';
    model.rhs = [ sum(xt); net.F(lim) + fx; -net.F(lim) + fx ];
    model.lb = [ zeros(1, 2*N + G), -inf(1, C) ]';
    model.ub = [ max(pdmax,0), max(pcmax,0), inf(1, G), zeros(1, C) ]';
    r = gurobi(model, params);
    assert(strcmp(r.status, 'OPTIMAL'), 'Step %d QP status %s', t, r.status);
    z = r.x';
    pd(:,t) = z(1:N); pc(:,t) = z(N+1:2*N);
    sg(:,t) = z(2*N+(1:G)); sc(:,t) = z(2*N+G+(1:C));
    ep = ep - pd(:,t)' ./ etad * dt + pc(:,t)' .* etac * dt;
    ep = min(max(ep, 0), par.E);
    e(:,t) = ep';
    lt = r.pi(1) * ones(1, net.nb);
    for kk = 1:nL, lt = lt + (r.pi(1+kk) + r.pi(1+nL+kk)) * net.PTDF(lim(kk), :); end
    lmp(t,:) = lt / dt;
    if useEma, Cm = rE .* Cm + (1 - rE) .* lmp(t, net.unitBus); end
end
tt = toc;
out.name = name; out.pd = pd; out.pc = pc; out.p = pd - pc; out.sg = sg; out.sc = sc; out.e = e;
out.lmp = lmp; out.v = vt; out.cost = net5_cost(pd, pc, sg, sc, par, net); out.ms = tt/T*1e3;
fprintf('[%-28s] cost %10.1f $ | %.2f ms/step\n', name, out.cost, out.ms);
end
