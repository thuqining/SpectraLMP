function res = netrun(cfg)
% NETRUN  One network experiment under the deployment protocol, for the
% cluster (one Slurm task per cfg). Fields of cfg (defaults in brackets):
%   tag       output name, saved to results/net/<tag>.mat
%   sys       '5bus' (PJM five-bus, net5_setup) or 'ieee14' (net14_setup)
%             ['5bus']
%   h         unit durations, h [1 4 30 200] (the base fleet uses
%             params_fleet exactly; other fleets follow fleet_refit)
%   F         limit of line E-D (5bus) or 5-6 (ieee14) in kW [8.6]; inf
%             for no limit
%   hist      initial history in days [60]
%   refit     refit period in days [30]
%   a, qs     cost of the generator at A (5bus) or of the units at buses
%             1, 2 (ieee14) [0.40 / 0.30] and q_+ [0.005]
%   scale     fleet power and energy multiplier [1]
%   eps       risk level rho [0.32]
%   methods   cell of 'prop','adp','mpc','myop','blind','cut','rl',
%             'hindavg','mpcexact','const','common','fixed'
%             [{'prop','adp'}]
% Every method is scored on the operating cost (1a) against the network
% hindsight optimum (1) of the same configuration.
gp = getenv('GUROBI_HOME'); cs = getenv('CPLEX_SHIM');       % cluster: CPLEX through cplexshim/gurobi.m
if ~isempty(cs), addpath(cs);
elseif ~isempty(gp), addpath(fullfile(gp, 'matlab')); else, addpath('C:\gurobi1302\win64\matlab'); end
here = fileparts(mfilename('fullpath')); root = fileparts(here); addpath(here);
if ~isfield(cfg, 'sys'), cfg.sys = '5bus'; end
is14 = strcmp(cfg.sys, 'ieee14');
d = struct('tag', 'test', 'h', [1 4 30 200], 'F', 8.6, 'hist', 60, 'refit', 30, ...
    'a', 0.40, 'qs', 0.005, 'scale', 1, 'eps', 0.32, 'trend', 7, 'eref', 0.6, 'termval', false, 'methods', {{'prop', 'adp'}});
if is14, d.a = 0.30; d.F = 3.1; end
fn = fieldnames(d);
for k = 1:numel(fn), if ~isfield(cfg, fn{k}), cfg.(fn{k}) = d.(fn{k}); end, end
spd = 96;
[~, dt, meta] = letter_data(fullfile(root, 'data', 'household_data_15min.csv'));
hrs = mod(hour(meta.t0) + minute(meta.t0)/60 + (0:size(meta.xb,1)-1)'*dt, 24);
if cfg.trend > 0
    Xs = meta.xb - movmean(meta.xb, [cfg.trend*spd-1, 0], 1);       % slow resource: causal trend-day mean
else
    Xs = meta.xb - mean(meta.xb, 1);                             % trend = 0: keep the seasons, remove only the mean
end
nTrain = 60*spd;
Ls = meta.xb;                                                  % net loads (for forecast errors)
if is14, [Xs, Ls] = net14_loads(Xs, net14_setup(), nTrain, meta.xb); end   % eleven load buses
xs = sum(Xs, 2); Xtest = Xs(nTrain+1:end, :); xtest = xs(nTrain+1:end); T = numel(xtest);
i0 = nTrain - cfg.hist*spd + 1; n0 = cfg.hist*spd;          % commissioning
Xsub = Xs(i0:end, :); xsub = xs(i0:end); hsub = hrs(i0:end);

% ---- fleet and costs
par = params_fleet(xs, dt);
h = cfg.h; N = numel(h);
if ~isequal(h, [1 4 30 200])
    hA = [1, 4, 30, 200]; c0A = [0.001, 0.003, 0.010, 0.030];
    etcA = [0.96, 0.95, 0.85, 0.65]; etdA = [0.96, 0.95, 0.85, 0.55];
    lh = @(v) log(min(max(v,1),200));
    w = 0.75.^(0:N-1); w = w/sum(w);
    par.N = N; par.h = h; par.P = 137.8*w; par.E = par.P.*h;
    par.c0 = exp(interp1(log(hA), log(c0A), lh(h))); par.etac = interp1(log(hA), etcA, lh(h));
    par.etad = interp1(log(hA), etdA, lh(h)); par.c = par.c0 + 0.01*(1 - par.etac.*par.etad);
    par.gamma = 0.1./h; par.name = arrayfun(@(v) sprintf('%gh', v), h, 'UniformOutput', false);
end
par.a_fuel = cfg.a; par.qs = cfg.qs; par.q = cfg.a/(2*std(xs(1:nTrain))) + cfg.qs;
par.eps_risk = cfg.eps; par.P = cfg.scale*par.P; par.E = cfg.scale*par.E;
par.eref = cfg.eref*par.E; par.e0 = par.eref;          % SoC reference and initial SoC
par.termval = cfg.termval;                             % score energy drawn from the initial stock

% ---- network: units spread over the candidate buses by duration rank
% (5bus: B, C, D, E; ieee14: 3, 4, 13, 2, i.e. the 30-h unit behind line 5-6)
if is14
    Fv = inf(20, 1); Fv(10) = cfg.F;                          % line 5-6
    net = net14_setup(Fv); net.genA(1:2) = cfg.a; cand = [3 4 13 2];
else
    Fv = inf(6, 1); Fv(6) = cfg.F;                            % line E-D
    net = net5_setup(Fv); net.genA(1) = cfg.a; cand = [2 3 4 5];
end
[~, rk] = sort(h); bus = zeros(1, N);
if N == 1, bus = cand(end); else, bus(rk) = cand(1 + round(3*(0:N-1)/(N-1))); end
net.unitBus = bus;

ev = evalc('pol = synth_policy(xsub(1:n0), par);'); %#ok<NASGU>
ev = evalc('Ho = method_hindsight_net(Xtest, par, net, struct(''terminalFree'', ~cfg.termval));'); %#ok<NASGU>
res = struct('cfg', cfg, 'unitBus', bus, 'lb', Ho.cost, 'congested', Ho.congested, ...
    'spread', mean(max(Ho.lmp,[],2) - min(Ho.lmp,[],2)));
if is14, res.flowP = prctile(abs(Ho.flow(10,:)), [50 65 80]); end   % |flow| of line 5-6 in the hindsight optimum
M = cfg.methods; has = @(m) any(strcmp(M, m));

% ---- monthly refit on all past data (network hindsight QP, 16-day hold)
needRefit = any(ismember(M, {'prop', 'adp', 'blind', 'const', 'common', 'hindcal', 'direct', 'adpf', 'split'}));
if has('hindcal')                       % hindsight averages of the test-year bus LMPs
    frh = (0:T-1)'/T; frh(frh > 0.5) = frh(frh > 0.5) - 1; wrh = 2*pi*frh/dt;
    Vh = zeros(T, N);
    for i = 1:N
        Vh(:,i) = real(ifft(fft(Ho.lmp(:, bus(i))) .* (pol.kappa(i)./(pol.kappa(i) + pol.ceff(i)*wrh.^2))));
    end
end
blk = cfg.refit*spd; nb = ceil(T/blk);
if needRefit
    V = struct('prop', zeros(T,N), 'blind', zeros(T,N), 'const', zeros(T,N), 'common', zeros(T,N), 'hindcal', zeros(T,N), ...
        'direct', zeros(T,N), 'adpfA', zeros(T,N));
    fits = cell(1, nb); adpfV1 = zeros(nb, N);
    VS = zeros(T, N); GS = zeros(nb, N); MS = zeros(nb, N); CS0 = zeros(1, N);
    tauU = sqrt(pol.ceff(:)'./pol.kappa(:)');
    if has('common'), pc = pol; tc = 34; pc.kappa = pol.ceff/tc^2; end   % one tau for all units
    for b = 1:nb
        b0 = (b-1)*blk + 1; b1 = min(T, b*blk); t0 = n0 + b0 - 1;
        Wn = Xsub(1:t0, :); w = xsub(1:t0);
        ev = evalc('Ht = method_hindsight_net(Wn, par, net);'); %#ok<NASGU>
        fits{b} = fit_values(w, par, Ht);
        o = struct('lossCal', true, 'Ht', Ht);
        if isfield(cfg, 'memH'), o.memH = cfg.memH; end
        if t0 < 16*spd, o.constant = true; end
        xf = xsub(t0+1:n0+b1); hf = hsub(1:n0+b1);
        pp = design_price(w, hsub(1:t0), par, pol, o);
        V.prop(b0:b1,:) = predict_values(xf, w, hf, pp);
        if has('hindcal')               % same calibration and fallback, exact future averages
            V.hindcal(b0:b1,:) = pp.cal(1,:) + pp.cal(2,:).*Vh(b0:b1,:);
        end
        if has('blind')
            net0 = net; net0.F(:) = inf;                       % same network, lines ignored
            ev = evalc('H1 = method_hindsight_net(Wn, par, net0);'); %#ok<NASGU>
            o1 = o; o1.Ht = H1;
            V.blind(b0:b1,:) = predict_values(xf, w, hf, design_price(w, hsub(1:t0), par, pol, o1));
        end
        if has('const')
            oc = o; oc.constant = true;
            V.const(b0:b1,:) = predict_values(xf, w, hf, design_price(w, hsub(1:t0), par, pol, oc));
        end
        if has('common')
            V.common(b0:b1,:) = predict_values(xf, w, hf, design_price(w, hsub(1:t0), par, pc, o));
        end
        if has('direct')                % ablation: same features and ridge, hindsight duals as targets
            od = o; od.direct = true;
            V.direct(b0:b1,:) = predict_values(xf, w, hf, design_price(w, hsub(1:t0), par, pol, od));
        end
        if has('adpf')                  % ADP with the same features plus SoC
            mH = min(192, floor(t0*dt/2)); if isfield(cfg, 'memH'), mH = min(cfg.memH, floor(t0*dt/2)); end
            pa = fit_adpf(w, hsub(1:t0), par, Ht, mH);
            Fa = price_features([w; xf], hf, mH); Fa = Fa(t0+1:end, :);
            V.adpfA(b0:b1,:) = [ones(size(Fa,1),1), (Fa - pa.mu)./pa.sd] * pa.B;
            adpfV1(b,:) = pa.v1;
        end
        if has('split')                 % causal LMP average observed + anticausal part predicted
            mH = min(192, floor(t0*dt/2)); if isfield(cfg, 'memH'), mH = min(cfg.memH, floor(t0*dt/2)); end
            ps = fit_split(w, hsub(1:t0), par, Ht, mH, tauU, bus);
            Fs = price_features([w; xf], hf, mH); Fs = Fs(t0+1:end, :);
            VS(b0:b1,:) = [ones(size(Fs,1),1), (Fs - ps.mu)./ps.sd] * ps.B;
            GS(b,:) = ps.g ./ ps.sdC; MS(b,:) = ps.muC;
            if b == 1, CS0 = ps.Cend; end
        end
    end
end
run = @(Vm, nm) method_values_net(Xtest, par, net, @(t, e) Vm(t,:), nm);
if has('prop'),   [~, R] = evalc('run(V.prop, ''prop'')');     res.prop = R.cost; res.msProp = R.ms; end
if has('prop') && ~has('split') && isfield(cfg, 'traj') && cfg.traj   % trajectories for Fig. 1(a)
    res.traj = struct('e', R.e, 'p', R.p, 'lmp', R.lmp, 'V', V.prop, 'Xtest', Xtest, 'E', par.E, 'P', par.P, ...
        'eH', Ho.e, 'muH', Ho.shadow_value);
end
if has('blind'),  [~, R] = evalc('run(V.blind, ''blind'')');   res.blind = R.cost; end
if has('direct'), [~, R] = evalc('run(V.direct, ''direct'')'); res.direct = R.cost; res.msDirect = R.ms; end
if has('split')
    [~, R] = evalc(['method_values_net(Xtest, par, net, @(t, e, Cm) VS(t,:) + GS(ceil(t/blk),:).*(Cm - MS(ceil(t/blk),:)), ' ...
        '''split'', struct(''emaTau'', tauU, ''emaInit'', CS0))']);
    res.split = R.cost; res.msSplit = R.ms;
    if isfield(cfg, 'traj') && cfg.traj   % trajectories of the proposed (split) dispatch for Figs. 1 and 2(c)
        res.traj = struct('e', R.e, 'p', R.p, 'lmp', R.lmp, 'V', R.v', 'Xtest', Xtest, 'E', par.E, 'P', par.P, ...
            'eH', Ho.e, 'muH', Ho.shadow_value);
    end
    res.lmpChk = [mean(Ho.lmp(:)), mean(R.lmp(:))];
end
if has('adpf')
    [~, R] = evalc('method_values_net(Xtest, par, net, @(t, e) V.adpfA(t,:) + adpfV1(ceil(t/blk),:).*e, ''adpf'')');
    res.adpf = R.cost; res.msAdpf = R.ms;
end
if has('const'),  [~, R] = evalc('run(V.const, ''const'')');   res.const = R.cost; end
if has('common'), [~, R] = evalc('run(V.common, ''common'')'); res.common = R.cost; end
if has('adp')
    [~, R] = evalc('method_values_net(Xtest, par, net, @(t, e) fits{ceil(t/blk)}(1,:) + fits{ceil(t/blk)}(2,:).*e, ''adp'')');
    res.adp = R.cost; res.msAdp = R.ms;
end
if has('fixed')                                               % design at commissioning, no refit
    ev = evalc('H0 = method_hindsight_net(Xsub(1:n0,:), par, net);'); %#ok<NASGU>
    o = struct('lossCal', true, 'Ht', H0); if n0 < 16*spd, o.constant = true; end
    pp = design_price(xsub(1:n0), hsub(1:n0), par, pol, o);
    Vf = predict_values(xtest, xsub(1:n0), hsub, pp);
    [~, R] = evalc('run(Vf, ''fixed'')'); res.fixed = R.cost;
end
if has('hindcal'), [~, R] = evalc('run(V.hindcal, ''hindcal'')'); res.hindcal = R.cost; end
if has('hindavg')                                             % hindsight averages of the bus LMPs
    fr = (0:T-1)'/T; fr(fr > 0.5) = fr(fr > 0.5) - 1; wr = 2*pi*fr/dt;
    Vo = zeros(T, N);
    for i = 1:N
        Vo(:,i) = real(ifft(fft(Ho.lmp(:, bus(i))) .* (pol.kappa(i)./(pol.kappa(i) + pol.ceff(i)*wr.^2))));
    end
    [~, R] = evalc('run(Vo, ''hindavg'')'); res.hindavg = R.cost;
end

% ---- benchmarks that set storage power directly
if any(ismember(M, {'myop', 'cut', 'rl', 'dqn'}))
    ev = evalc('polf = finish_policy(pol, xsub(1:n0), par, 1);'); %#ok<NASGU>
end
if has('myop')
    g = 1 ./ par.c; wg = g / (1/par.q + sum(g));
    [~, R] = evalc('method_policy_net(Xtest, par, net, @(t, e) wg*(xtest(t) - polf.mu) + polf.gamma.*(e - par.eref), ''myop'')');
    res.myop = R.cost; res.msMyop = R.ms;
end
if has('cut')
    xfit = xsub(1:n0-25*spd); xval = xsub(n0-25*spd+1:n0);
    ev = evalc('p0 = synth_policy(xfit, par);'); %#ok<NASGU>
    p0.valueFit = fit_values(xfit, par);
    ev = evalc('wc = tune_cutoffs(xval, par, p0);'); %#ok<NASGU>
    ev = evalc('H60 = method_hindsight_net(Xsub(1:n0,:), par, net);'); %#ok<NASGU>
    f60 = fit_values(xsub(1:n0), par, H60);
    hist = polf.history(:)'; xz = [hist, xtest(:)'] - polf.mu;
    ma = @(v, W) filter(ones(1,W)/W, 1, v);
    m1 = ma(xz, wc(1)); m2 = ma(xz, wc(2)); m3 = ma(xz, wc(3));
    bd = [xz - m1; m1 - m2; m2 - m3; m3]; bd = bd(:, numel(hist)+1:end);
    [~, R] = evalc(['method_policy_net(Xtest, par, net, @(t, e) bd(:,t)'' + polf.gamma.*(e - par.eref), ''cut'', ' ...
        '@(t, e) f60(1,:) + f60(2,:).*e)']);
    res.cut = R.cost; res.msCut = R.ms; res.cutWin = wc;
end
if has('rl')
    Wrl = rl_train(xsub(1:n0), par, polf);
    sx = std(xsub(1:n0)); xzt = xtest(:)' - polf.mu;
    m4 = filter(1/4, [1, -0.75], xzt); m96 = filter(1/96, [1, -(1-1/96)], xzt);
    tod = 2*pi*mod((1:T)-1, 96)/96;
    rls = @(t, e) par.P .* tanh([1, xzt(t)/sx, m4(t)/sx, m96(t)/sx, (e./par.E - 0.5), sin(tod(t)), cos(tod(t))] * Wrl);
    [~, R] = evalc('method_policy_net(Xtest, par, net, rls, ''rl'')');
    res.rl = R.cost; res.msRl = R.ms;
end
if has('none')                                                % no storage dispatch (reference)
    [~, R] = evalc('method_policy_net(Xtest, par, net, @(t, e) zeros(1, N), ''none'')');
    res.none = R.cost;
end
if has('dqnm')                    % DQN [5] retrained at each monthly refit on all past data (dqn_train.py per block)
    dd = fullfile(root, 'results', 'dqnm'); if isfield(cfg, 'dqnmDir'), dd = cfg.dqnmDir; end
    if ~isfolder(dd), mkdir(dd); end
    P = par.P; E = par.E; etac = par.etac; etad = par.etad; c0 = par.c0; %#ok<NASGU>
    a_fuel = par.a_fuel; qs = par.qs; qc = par.qc; e0 = par.e0; %#ok<NASGU>
    MU = zeros(1, nb); SX = zeros(1, nb);
    for b = 1:nb
        t0 = n0 + (b-1)*blk; x = xsub(1:t0); mu = mean(x); sx = std(x); MU(b) = mu; SX(b) = sx;
        save(fullfile(dd, sprintf('env_%02d.mat', b)), 'x', 'P', 'E', 'etac', 'etad', 'c0', 'a_fuel', 'qs', 'qc', 'e0', 'mu', 'sx', 'dt');
    end
    if isfield(cfg, 'dqnEnvOnly') && cfg.dqnEnvOnly, disp('dqnm envs written'); return; end
    Qb = cell(1, nb); EXb = zeros(T, 5, nb); tod = 2*pi*mod((1:T)-1, 96)/96;
    for b = 1:nb
        S = load(fullfile(dd, sprintf('w_%02d.mat', b))); Q = struct('LEV', S.LEV(:)');
        for i = 1:N
            Q.W1{i} = S.(sprintf('W1_%d', i)); Q.b1{i} = S.(sprintf('b1_%d', i))(:)';
            Q.W2{i} = S.(sprintf('W2_%d', i)); Q.b2{i} = S.(sprintf('b2_%d', i))(:)';
            Q.W3{i} = S.(sprintf('W3_%d', i)); Q.b3{i} = S.(sprintf('b3_%d', i))(:)';
        end
        Qb{b} = Q; xzt = xtest(:)' - MU(b);
        m4 = filter(1/4, [1, -0.75], xzt); m96 = filter(1/96, [1, -(1-1/96)], xzt);
        EXb(:,:,b) = [min(max([xzt'/SX(b), m4'/SX(b), m96'/SX(b)], -5), 5), sin(tod'), cos(tod')];
    end
    [~, R] = evalc(['method_policy_net(Xtest, par, net, @(t, e) dqn_policy(Qb{ceil(t/blk)}, ' ...
        'EXb(t,:,ceil(t/blk)), e, par), ''dqnm'')']);
    res.dqnm = R.cost; res.msDqnm = R.ms;
end
if has('dqn')                                                 % double deep Q-learning [5], one agent per unit
    envf = [tempname '.mat']; outf = [tempname '.mat'];
    x = xsub(1:n0); P = par.P; E = par.E; etac = par.etac; etad = par.etad; c0 = par.c0; %#ok<NASGU>
    a_fuel = par.a_fuel; qs = par.qs; qc = par.qc; e0 = par.e0; mu = polf.mu; sx = std(xsub(1:n0)); %#ok<NASGU>
    save(envf, 'x', 'P', 'E', 'etac', 'etad', 'c0', 'a_fuel', 'qs', 'qc', 'e0', 'mu', 'sx', 'dt');
    if isfield(cfg, 'dqnW') && isfile(cfg.dqnW)               % reuse trained weights (same training data)
        outf = cfg.dqnW;
    else
        py = getenv('PYEXE'); if isempty(py), py = 'python'; end
        [st, msg] = system(sprintf('"%s" "%s" "%s" "%s"', py, fullfile(here, 'dqn_train.py'), envf, outf));
        disp(msg); assert(st == 0, 'dqn_train failed');
        copyfile(outf, fullfile(root, 'results', 'dqn_weights.mat'));
    end
    S = load(outf); Q = struct('LEV', S.LEV(:)');
    for i = 1:N
        Q.W1{i} = S.(sprintf('W1_%d', i)); Q.b1{i} = S.(sprintf('b1_%d', i))(:)';
        Q.W2{i} = S.(sprintf('W2_%d', i)); Q.b2{i} = S.(sprintf('b2_%d', i))(:)';
        Q.W3{i} = S.(sprintf('W3_%d', i)); Q.b3{i} = S.(sprintf('b3_%d', i))(:)';
    end
    xzt = xtest(:)' - polf.mu; sxt = std(xsub(1:n0));
    m4 = filter(1/4, [1, -0.75], xzt); m96 = filter(1/96, [1, -(1-1/96)], xzt);
    tod = 2*pi*mod((1:T)-1, 96)/96;
    EX = [min(max([xzt'/sxt, m4'/sxt, m96'/sxt], -5), 5), sin(tod'), cos(tod')];
    [~, R] = evalc('method_policy_net(Xtest, par, net, @(t, e) dqn_policy(Q, EX(t,:), e, par), ''dqn'')');
    res.dqn = R.cost; res.msDqn = R.ms; res.dqnGamma = S.gamma;
end
if has('mpc')
    [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''persistence'')'); res.mpc = R.cost; res.msMpc = R.ms;
end
sgm = std(Xsub(1:n0, :));                                     % training deviation by bus
for r = [10 20]                                               % uniform +-r% of the net load
    m = sprintf('mpcu%d', r);
    if has(m)
        Lt = Ls(nTrain+1:end, :);
        [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''uniform'', struct(''err'', r/100, ''sigma'', 0, ''L'', Lt))');
        res.(m) = R.cost;
    end
end
for r = [5 10 20 50]
    m = sprintf('mpc%d', r);
    if has(m)
        [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''noisy'', struct(''err'', r/100, ''sigma'', sgm))');
        res.(m) = R.cost;
    end
end
if has('mpcfit')                                              % data-driven forecast, past data only
    Xt = Xsub(1:n0, :); nl = size(Xt, 2); Hh = 96; coef = zeros(3, Hh, nl); rmse = zeros(Hh, nl);
    for hh = 1:Hh
        tt = (97 - hh):(n0 - hh); tt = tt(tt >= 1);
        for k = 1:nl
            A = [Xt(tt, k), Xt(tt + hh - 96, k), ones(numel(tt), 1)]; y = Xt(tt + hh, k);
            coef(:, hh, k) = A \ y; rmse(hh, k) = std(y - A*coef(:, hh, k));
        end
    end
    w = sgm > 0; res.fitErr1 = mean(rmse(1, w) ./ sgm(w)); res.fitErr96 = mean(rmse(96, w) ./ sgm(w));
    Xh0 = Xsub(n0-95:n0, :);                                  % last training day
    [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''fit'', struct(''coef'', coef, ''Xhist'', Xh0))');
    res.mpcfit = R.cost;
end
if has('mpcxgb')                                              % XGBoost forecast from past data only
    Fx = load(fullfile(root, 'results', 'forecast_xgb_s.mat'));
    assert(size(Fx.Fg, 1) == T && size(Fx.Map, 2) == size(Xtest, 2) && cfg.hist == 60 && cfg.trend == 0);
    [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''given'', struct(''Fg'', Fx.Fg, ''Map'', Fx.Map))');
    res.mpcxgb = R.cost; res.msMpcxgb = R.ms;
end
if any(startsWith(M, 'mpc'))                                  % 24-h error of the persistence forecast
    w = sgm > 0;
    res.persErr = mean(std(Xtest(97:end, w) - Xtest(1:end-96, w)) ./ sgm(w));
end
if has('mpcexact')
    [~, R] = evalc('method_mpc_net(Xtest, par, net, pol, ''perfect'')'); res.mpcexact = R.cost;
end
od = fullfile(root, 'results', 'net'); if ~isfolder(od), mkdir(od); end
save(fullfile(od, [cfg.tag '.mat']), 'res');
f = setdiff(fieldnames(res), {'cfg', 'unitBus', 'lb', 'congested', 'spread', 'cutWin', 'traj', 'flowP', 'lmpChk', 'dqnGamma'});
if isfield(res, 'lmpChk'), fprintf('mean LMP hindsight %.4f online %.4f\n', res.lmpChk); end
s = sprintf('[%s] F %g congested %.1f%% |', cfg.tag, cfg.F, 100*res.congested);
for k = 1:numel(f)
    if ~startsWith(f{k}, 'ms'), s = [s sprintf(' %s %.1f', f{k}, 100*res.(f{k})/res.lb)]; end %#ok<AGROW>
end
disp(s);
end

function W = rl_train(xtrain, par, pol)
% CEM policy search as method_rl_cem (seed 7, 40 x 30).
N = par.N; sx = std(xtrain); nf = 4 + N + 2;
rng(7); npop = 40; nelite = 8; niter = 30;
Wm = zeros(N*nf, 1); Ws = 0.3*ones(N*nf, 1); xtr = xtrain(:)';
for it = 1:niter
    cand = Wm + Ws .* randn(N*nf, npop); Jc = zeros(1, npop);
    for k = 1:npop, Jc(k) = rollout(reshape(cand(:,k), nf, N), xtr, par, pol, sx); end
    [~, ord] = sort(Jc); el = cand(:, ord(1:nelite));
    Wm = mean(el, 2); Ws = std(el, 0, 2) + 0.02;
end
W = reshape(Wm, nf, N);
end

function J = rollout(W, x, par, pol, sx)
T = numel(x); N = par.N; p = zeros(N,T); s = zeros(1,T);
e_prev = par.e0(:)'; ma4 = 0; ma96 = 0;
for t = 1:T
    xz = x(t) - pol.mu; ma4 = 0.75*ma4 + 0.25*xz; ma96 = (1-1/96)*ma96 + xz/96;
    tod = 2*pi*mod(t-1, 96)/96;
    phi = [1; xz/sx; ma4/sx; ma96/sx; (e_prev(:)./par.E(:) - 0.5); sin(tod); cos(tod)];
    [pt, st, e_prev] = util_clip(par.P .* tanh(phi' * W), x(t), e_prev, par);
    p(:,t) = pt'; s(t) = st;
end
J = sum(par.a_fuel*max(s,0) + par.qs*max(s,0).^2 + par.qc*min(s,0).^2 + par.c0*p.^2) * par.dt;
end
