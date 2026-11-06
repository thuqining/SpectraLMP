function out = method_policy_net(X, par, net, stepFun, name, valueFun)
% METHOD_POLICY_NET  Network evaluation of a benchmark that sets storage
% power directly (filters, RL, myopic rule). At each step the benchmark
% requests p0 = stepFun(t, e) from the system net load, clipped to the
% power and energy limits, and a one-step DC-OPF QP dispatches the rest:
%   no valueFun : storage follows p0, with a large tracking penalty so
%                 that it deviates only where the line limits require;
%   valueFun    : completion as in the single-bus case -- each unit may
%                 move from p0 only toward the residual sum(x)-sum(p0),
%                 priced by its stored-energy value v = valueFun(t, e).
%   X  T x numel(net.loadBus) nodal net loads (test period)
if nargin < 6, valueFun = []; end
[T, ~] = size(X); N = par.N; dt = par.dt;
G = numel(net.genBus); C = numel(net.curtBus);
lim = find(isfinite(net.F)); nL = numel(lim);
Gu = net.PTDF(lim, net.unitBus); Gg = net.PTDF(lim, net.genBus);
Gc = net.PTDF(lim, net.curtBus); Gx = net.PTDF(lim, net.loadBus);
c0 = par.c0(:)'; etac = par.etac(:)'; etad = par.etad(:)';
Wt = 10;                                       % $/(kW^2 h) tracking weight
Arow = [ones(1,N), -ones(1,N), ones(1,G), ones(1,C)];
Al = [Gu, -Gu, Gg, Gc];
Ad = [eye(N), -eye(N), zeros(N, G+C)];         % net power p = pd - pc
params.OutputFlag = 0; params.Threads = 1;
if isempty(valueFun)
    Q = blkdiag([diag(c0+Wt), diag(c0-Wt); diag(c0-Wt), diag(c0+Wt)], diag(net.genQ), net.qc*eye(C));
    model.A = sparse([Arow; Al; Al]);
    model.sense = [ '='; repmat('<', nL, 1); repmat('>', nL, 1) ];
else
    Q = blkdiag([diag(c0), diag(c0); diag(c0), diag(c0)], diag(net.genQ), net.qc*eye(C));
end
model.Q = sparse(Q) * dt;
pd = zeros(N,T); pc = zeros(N,T); sg = zeros(G,T); sc = zeros(C,T); e = zeros(N,T);
ep = par.e0(:)';
tic;
for t = 1:T
    xt = X(t,:);
    pdmax = min(par.P, etad .* ep / dt);
    pcmax = min(par.P, (par.E - ep) ./ (etac * dt));
    p0 = min(pdmax, max(-pcmax, stepFun(t, ep)));
    fx = Gx * xt';
    rhs = [ sum(xt); net.F(lim) + fx; -net.F(lim) + fx ];
    model.lb = [ zeros(1, 2*N + G), -inf(1, C) ]';
    model.ub = [ max(pdmax,0), max(pcmax,0), inf(1, G), zeros(1, C) ]';
    if isempty(valueFun)
        model.obj = [ -2*Wt*p0*dt, 2*Wt*p0*dt, net.genA(:)'*dt, zeros(1, C) ]';
        model.rhs = rhs;
    else
        v = max(valueFun(t, ep), 0);
        model.obj = [ (v ./ etad) * dt, -(v .* etac) * dt, net.genA(:)' * dt, zeros(1, C) ]';
        delta = sum(xt) - sum(p0);
        if delta >= 0, sd = '>'; else, sd = '<'; end
        model.A = sparse([Arow; Al; Al; Ad]);
        model.sense = [ '='; repmat('<', nL, 1); repmat('>', nL, 1); repmat(sd, N, 1) ];
        model.rhs = [ rhs; p0(:) ];
    end
    r = gurobi(model, params);
    if ~any(strcmp(r.status, {'OPTIMAL', 'SUBOPTIMAL'}))
        error('Step %d QP status %s', t, r.status);
    end
    z = r.x';
    pd(:,t) = z(1:N); pc(:,t) = z(N+1:2*N);
    sg(:,t) = z(2*N+(1:G)); sc(:,t) = z(2*N+G+(1:C));
    ep = ep - pd(:,t)' ./ etad * dt + pc(:,t)' .* etac * dt;
    ep = min(max(ep, 0), par.E);
    e(:,t) = ep';
end
tt = toc;
out.name = name; out.pd = pd; out.pc = pc; out.p = pd - pc; out.sg = sg; out.sc = sc; out.e = e;
out.cost = net5_cost(pd, pc, sg, sc, par, net); out.ms = tt/T*1e3;
fprintf('[%-28s] cost %10.1f $ | %.2f ms/step\n', name, out.cost, out.ms);
end
