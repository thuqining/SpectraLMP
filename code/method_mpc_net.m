function out = method_mpc_net(X, par, net, pol, fcType, fo)
% METHOD_MPC_NET  Network MPC benchmark: 24-h horizon, 30-min recourse,
% a DC-OPF QP over the horizon with the soft terminal penalty
% kappa_i (e_H - e_ref)^2 of the single-bus MPC. fcType = 'persistence'
% (each bus repeats its net load of the previous day), 'perfect', or
% 'noisy': the true net load plus an error that grows with the lead time as
% a random walk, with standard deviation fo.err*fo.sigma(n) at 24 h
% (fo.sigma: training deviation of each bus; seed fo.seed, default 1), or
% 'fit': a data-driven forecast from past data only, x_hat(t0-1+h) =
% a_h x(t0-1) + b_h x(t0-1+h-96) + c_h per bus, with fo.coef (3 x H x nl)
% fitted on the training data (netrun, method 'mpcfit'), or
% 'uniform': the true value of each future step plus fo.err*u*fo.L, with u
% uniform on [-1, 1], independent across steps, buses, and solves, and
% fo.L (T x nl) the net load of each bus (seed fo.seed, default 1).
% The committed storage power is applied to the actual net load by a
% one-step DC-OPF in which generators, peakers, and curtailment adjust.
[T, nl] = size(X); N = par.N; dt = par.dt; H = 96; rc = 2;
G = numel(net.genBus); C = numel(net.curtBus);
lim = find(isfinite(net.F)); nL = numel(lim);
Gu = net.PTDF(lim, net.unitBus); Gg = net.PTDF(lim, net.genBus);
Gc = net.PTDF(lim, net.curtBus); Gx = net.PTDF(lim, net.loadBus);
etac = par.etac(:); etad = par.etad(:); c0 = par.c0(:);
NH = N*H; GH = G*H; CH = C*H; Ih = speye(H);
o_pd = 0; o_pc = NH; o_e = 2*NH; o_sg = 3*NH; o_sc = 3*NH + GH; nv = 3*NH + GH + CH;
% ---- horizon QP structure (rhs and bounds change per solve) ----
Qpd = kron(Ih, sparse(diag(c0))) * dt;
Q = blkdiag(Qpd, Qpd, sparse(NH, NH), kron(Ih, sparse(diag(net.genQ))) * dt, ...
            kron(Ih, sparse(net.qc*eye(C))) * dt);
Q(o_pd+(1:NH), o_pc+(1:NH)) = Qpd; Q(o_pc+(1:NH), o_pd+(1:NH)) = Qpd;
termI = o_e + (H-1)*N + (1:N);
kap = pol.kappa(:);
Q = Q + sparse(termI, termI, kap, nv, nv);
obj = zeros(nv, 1); obj(o_sg+(1:GH)) = repmat(net.genA(:), H, 1) * dt;
obj(termI) = -2 * kap .* par.eref(:);
Ab = [kron(Ih, sparse(ones(1,N))), -kron(Ih, sparse(ones(1,N))), sparse(H, NH), ...
      kron(Ih, sparse(ones(1,G))), kron(Ih, sparse(ones(1,C)))];
As = [kron(Ih, sparse(diag(dt./etad))), -kron(Ih, sparse(diag(etac*dt))), ...
      speye(NH) - spdiags(ones(NH,1), -N, NH, NH), sparse(NH, GH + CH)];
Al = sparse(0, nv);
for kk = 1:nL
    Al = [Al; kron(Ih, sparse(Gu(kk,:))), -kron(Ih, sparse(Gu(kk,:))), sparse(H, NH), ...
          kron(Ih, sparse(Gg(kk,:))), kron(Ih, sparse(Gc(kk,:)))]; %#ok<AGROW>
end
M.Q = Q; M.obj = obj; M.A = [Ab; As; Al; Al];
M.sense = [repmat('=', H + NH, 1); repmat('<', nL*H, 1); repmat('>', nL*H, 1)];
M.lb = [zeros(3*NH + GH, 1); -inf(CH, 1)];
M.ub = [repmat(par.P(:), H, 1); repmat(par.P(:), H, 1); repmat(par.E(:), H, 1); inf(GH, 1); zeros(CH, 1)];
prm.OutputFlag = 0; prm.Threads = 4;
% ---- one-step dispatch with storage fixed ----
S.Q = sparse(blkdiag(diag(net.genQ), net.qc*eye(C))) * dt;
S.obj = [net.genA(:) * dt; zeros(C, 1)];
S.A = sparse([ones(1, G + C); [Gg, Gc]; [Gg, Gc]]);
S.sense = ['='; repmat('<', nL, 1); repmat('>', nL, 1)];
S.lb = [zeros(G,1); -inf(C,1)]; S.ub = [inf(G,1); zeros(C,1)];
pd = zeros(N,T); pc = zeros(N,T); sg = zeros(G,T); sc = zeros(C,T); e = zeros(N,T);
ep = par.e0(:); tsum = 0; ns = 0; xmean = mean(X, 1);
useHist = nargin > 5 && isfield(fo, 'Xhist');         % past net load before the test (no test data)
if useHist, XA = [fo.Xhist; X]; oh = size(fo.Xhist, 1); end
if any(strcmp(fcType, {'noisy', 'uniform'}))
    if ~isfield(fo, 'seed'), fo.seed = 1; end
    rs0 = RandStream('mt19937ar', 'Seed', fo.seed);
    sw = fo.err * fo.sigma(:)' / sqrt(H);           % random-walk step, 1 x nl
end
for t0 = 1:rc:T
    idx = t0:t0+H-1;
    switch fcType
        case 'perfect',     Xh = X(min(idx, T), :);
        case 'noisy',       Xh = X(min(idx, T), :) + cumsum(randn(rs0, H, nl) .* sw, 1);
        case 'uniform',     Xh = X(min(idx, T), :) + fo.err * (2*rand(rs0, H, nl) - 1) .* fo.L(min(idx, T), :);
        case 'fit'
            if useHist                              % lags from the actual past, never the test mean
                xl = XA(oh + t0 - 1, :); Xy = XA(oh + t0 - 1 + (1:H)' - 96, :);
            else
                if t0 > 1, xl = X(t0-1, :); else, xl = xmean; end
                iy = t0 - 1 + (1:H)' - 96; Xy = repmat(xmean, H, 1); m = iy >= 1; Xy(m,:) = X(iy(m), :);
            end
            Xh = squeeze(fo.coef(1,:,:)) .* xl + squeeze(fo.coef(2,:,:)) .* Xy + squeeze(fo.coef(3,:,:));
        case 'given'                                % precomputed causal forecast (forecast_xgb.py)
            Xh = double(reshape(fo.Fg(t0, :, :), H, [])) * fo.Map;
        case 'persistence', ip = idx - 96; Xh = repmat(xmean, H, 1);
                            m = ip >= 1; Xh(m,:) = X(ip(m), :);
    end
    rs = zeros(NH, 1); rs(1:N) = ep;
    fx = Xh * Gx';                                  % H x nL
    ru = reshape(net.F(lim)' + fx, [], 1);
    rl = reshape(-net.F(lim)' + fx, [], 1);
    M.rhs = [sum(Xh, 2); rs; ru; rl];
    tic; r = gurobi(M, prm); tsum = tsum + toc; ns = ns + 1;
    if ~any(strcmp(r.status, {'OPTIMAL', 'SUBOPTIMAL'})) || ~isfield(r, 'x')
        error('MPC QP status %s at %d', r.status, t0);
    end
    z = r.x;
    for j = 1:rc
        t = t0 + j - 1; if t > T, break; end
        qd = z(o_pd + (j-1)*N + (1:N)); qc = z(o_pc + (j-1)*N + (1:N));
        qd = min(qd, etad .* ep / dt); qc = min(qc, (par.E(:) - ep) ./ (etac * dt));
        qd = max(qd, 0); qc = max(qc, 0);
        pnet = qd - qc;
        xt = X(t, :);
        S.rhs = [sum(xt) - sum(pnet); net.F(lim) + Gx*xt' - Gu*pnet; -net.F(lim) + Gx*xt' - Gu*pnet];
        rr = gurobi(S, prm);
        if ~strcmp(rr.status, 'OPTIMAL')            % storage set point infeasible: drop it
            qd(:) = 0; qc(:) = 0; pnet = zeros(N,1);
            S.rhs = [sum(xt); net.F(lim) + Gx*xt'; -net.F(lim) + Gx*xt'];
            rr = gurobi(S, prm);
        end
        pd(:,t) = qd; pc(:,t) = qc; sg(:,t) = rr.x(1:G); sc(:,t) = rr.x(G+1:end);
        ep = min(max(ep - qd ./ etad * dt + qc .* etac * dt, 0), par.E(:));
        e(:,t) = ep;
    end
end
out.name = sprintf('MPC (%s, 24 h, network)', fcType);
out.pd = pd; out.pc = pc; out.p = pd - pc; out.sg = sg; out.sc = sc; out.e = e;
out.cost = net5_cost(pd, pc, sg, sc, par, net); out.ms = tsum/ns*1e3/rc;
fprintf('[%-28s] cost %10.1f $ | %.2f ms/step\n', out.name, out.cost, out.ms);
end
