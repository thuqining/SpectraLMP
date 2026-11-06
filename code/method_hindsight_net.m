function out = method_hindsight_net(X, par, net, opts)
% METHOD_HINDSIGHT_NET  Network version of problem (1): full-horizon
% hindsight DC-OPF with storage, generators, curtailment, and PTDF line
% limits. Simultaneous charging and discharging is allowed (relaxation).
%   X     T x numel(net.loadBus) net loads at the load buses (kW)
%   opts.terminalFree  true: no terminal value (oracle, operating cost);
%                      false (default): energy left valued at a*etad
%                      (training QP, as method_hindsight)
% Outputs: dispatch, operating cost (1a), balance price, nodal LMPs, the
% LMP at each unit's bus (priceUnit), SoC duals, and line flows.
% Stacking z = [pd(NT); pc(NT); e(NT); sg(GT); sc(CT)], (i,t) -> (t-1)N+i.
if nargin < 4, opts = struct(); end
free = isfield(opts, 'terminalFree') && opts.terminalFree;
[T, nl] = size(X); N = par.N; dt = par.dt;
G = numel(net.genBus); C = numel(net.curtBus);
NT = N*T; GT = G*T; CT = C*T;
o_pd = 0; o_pc = NT; o_e = 2*NT; o_sg = 3*NT; o_sc = 3*NT + GT;
nv = 3*NT + GT + CT;
It = speye(T);

% ---------- objective ----------
c0 = par.c0(:)';                                        % c0 (pd+pc)^2 per step
Qpd = kron(It, sparse(diag(c0))) * dt;
Q = blkdiag(Qpd, Qpd, sparse(NT, NT), kron(It, sparse(diag(net.genQ))) * dt, ...
            kron(It, sparse(net.qc*eye(C))) * dt);
Q(o_pd+(1:NT), o_pc+(1:NT)) = Qpd;                      % cross term 2 c0 pd pc (split)
Q(o_pc+(1:NT), o_pd+(1:NT)) = Qpd;
obj = zeros(nv, 1);
obj(o_sg+(1:GT)) = repmat(net.genA(:), T, 1) * dt;
terminalI = o_e + (T-1)*N + (1:N);
if ~free, obj(terminalI) = -par.a_fuel * par.etad(:); end

% ---------- balance ----------
Ab = [kron(It, sparse(ones(1,N))), -kron(It, sparse(ones(1,N))), sparse(T, NT), ...
      kron(It, sparse(ones(1,G))), kron(It, sparse(ones(1,C)))];
rb = sum(X, 2);

% ---------- SoC recursion ----------
etac = par.etac(:); etad = par.etad(:);
As = [kron(It, sparse(diag(dt./etad))), -kron(It, sparse(diag(etac*dt))), ...
      speye(NT) - spdiags(ones(NT,1), -N, NT, NT), sparse(NT, GT + CT)];
rs = zeros(NT, 1); rs(1:N) = par.e0(:);

% ---------- line limits ----------
lim = find(isfinite(net.F));
Gu = net.PTDF(:, net.unitBus); Gg = net.PTDF(:, net.genBus);
Gc = net.PTDF(:, net.curtBus); Gx = net.PTDF(:, net.loadBus);
Al = sparse(0, nv); ru = []; rl = [];
for l = lim(:)'
    row = [kron(It, sparse(Gu(l,:))), -kron(It, sparse(Gu(l,:))), sparse(T, NT), ...
           kron(It, sparse(Gg(l,:))), kron(It, sparse(Gc(l,:)))];
    Al = [Al; row]; %#ok<AGROW>
    ru = [ru; net.F(l) + X*Gx(l,:)']; %#ok<AGROW>
end
nL = numel(lim);

model.Q = Q; model.obj = obj;
model.A = [Ab; As; Al; Al];
model.rhs = [rb; rs; ru; ru - 2*kron(net.F(lim), ones(T,1))];
model.sense = [repmat('=', T + NT, 1); repmat('<', nL*T, 1); repmat('>', nL*T, 1)];
model.lb = [zeros(3*NT + GT, 1); -inf(CT, 1)];
model.ub = [repmat(par.P(:), T, 1); repmat(par.P(:), T, 1); repmat(par.E(:), T, 1); ...
            inf(GT, 1); zeros(CT, 1)];
params.OutputFlag = 0; params.Threads = 4;
tic; r = gurobi(model, params); tsol = toc;
assert(strcmp(r.status, 'OPTIMAL'), 'Network hindsight QP status: %s', r.status);
z = r.x;
pd = reshape(z(o_pd+(1:NT)), N, T); pc = reshape(z(o_pc+(1:NT)), N, T);
e  = reshape(z(o_e +(1:NT)), N, T);
sg = reshape(z(o_sg+(1:GT)), G, T); sc = reshape(z(o_sc+(1:CT)), C, T);

pib = r.pi(1:T);
lmp = repmat(pib, 1, net.nb);
if nL > 0
    k0 = T + NT;
    for kk = 1:nL
        l = lim(kk);
        piu = r.pi(k0 + (kk-1)*T + (1:T));
        pil = r.pi(k0 + nL*T + (kk-1)*T + (1:T));
        lmp = lmp + (piu + pil) * net.PTDF(l, :);
    end
end
out.pd = pd; out.pc = pc; out.p = pd - pc; out.e = e; out.sg = sg; out.sc = sc;
out.price = pib' / dt;                       % balance price (slack bus), $/kWh
out.lmp = lmp / dt;                          % T x nb nodal LMPs, $/kWh
out.priceUnit = out.lmp(:, net.unitBus);     % LMP at each unit's bus
out.shadow_value = -reshape(r.pi(T+(1:NT)), N, T);
inj = zeros(net.nb, T);
for i = 1:N, inj(net.unitBus(i),:) = inj(net.unitBus(i),:) + out.p(i,:); end
for g = 1:G, inj(net.genBus(g),:) = inj(net.genBus(g),:) + sg(g,:); end
for c = 1:C, inj(net.curtBus(c),:) = inj(net.curtBus(c),:) + sc(c,:); end
for n = 1:nl, inj(net.loadBus(n),:) = inj(net.loadBus(n),:) - X(:,n)'; end
out.flow = net.PTDF * inj;                   % L x T, kW
out.cost = net5_cost(pd, pc, sg, sc, par, net);
out.congested = mean(any(abs(out.flow(lim,:)) >= net.F(lim) - 1e-3, 1));
out.solveTime = tsol;
fprintf('[hindsight-net] T = %d, %.1f s, cost %.1f $, congested %.1f%% of steps\n', ...
    T, tsol, out.cost, 100*out.congested);
end
