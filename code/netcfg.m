function cfg = netcfg(id, set)
% NETCFG  Configurations of the five-bus study (one per Slurm array task).
%   set 'check'  reproduce net5_main at 8.6 kW (112.5 / 115.4 / 118.9)
%   set 'main'   Table I and the hindsight decomposition, F = inf and 8.6
%   set 'fleet'  six fleets x {inf, 8.6 kW} x {60, 30, 7} days of history
%   set 'sens'   twelve one-at-a-time settings x {inf, 8.6 kW}
%   set 'lim'    four methods at the 6.0 and 4.0 kW limits (Fig. 2)
% A set name ending in '14' (e.g. 'main14') runs the same design on the
% IEEE 14-bus system: line 5-6 at 3.1 kW (65th percentile of its
% unconstrained flow; 4.2 and 2.4 kW for the 80th and 50th), base a = 0.30.
% A set name ending in '14s' (e.g. 'main14s') is the seasonal version of the
% 14-bus study: the net load keeps its seasons (only its mean is removed),
% energy drawn from the initial stock is scored, and line 5-6 is limited to
% 3.3 kW (65th percentile of its unconstrained flow; 4.3 and 2.8 kW for the
% 80th and 50th). Tags start with 'i14s_'.
isS = endsWith(set, '14s');
if isS, set = [set(1:end-1)]; end
is14 = endsWith(set, '14');
if is14, set = set(1:end-2); end
FL = [inf, 8.6]; ft = {'free', 'lim'}; aL = [0.30 0.50]; F2 = [6.0, 4.0];
if is14, FL = [inf, 3.1]; aL = [0.25 0.35]; F2 = [4.2, 2.4]; end
if isS, FL = [inf, 3.3]; F2 = [4.3, 2.8]; end
switch set
    case 'check'
        cfg = struct('tag', 'check_q80', 'F', 8.6, 'methods', {{'prop', 'adp', 'myop'}});
    case 'main'
        k = id;                                       % 1..2
        cfg = struct('tag', ['main_' ft{k}], 'F', FL(k), 'methods', ...
            {{'prop', 'adp', 'mpc', 'myop', 'cut', 'rl', 'hindavg', 'mpcexact', ...
              'blind', 'const', 'common', 'fixed'}});
    case 'fleet'
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        hd = [60 30 7];
        [f, l, d] = ind2sub([6, 2, 3], id);           % 1..36
        m = {'prop', 'adp'}; if hd(d) == 60, m = {'prop', 'adp', 'mpc', 'myop'}; end
        cfg = struct('tag', sprintf('fleet_%s_%d_%d', ft{l}, hd(d), f), 'h', H{f}, 'F', FL(l), ...
            'hist', hd(d), 'methods', {m});
    case 'abl'                                        % direct-value ablation and feature ADP, same 36 settings as 'fleet'
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        hd = [60 30 7];
        [f, l, d] = ind2sub([6, 2, 3], id);           % 1..36
        cfg = struct('tag', sprintf('abl_%s_%d_%d', ft{l}, hd(d), f), 'h', H{f}, 'F', FL(l), ...
            'hist', hd(d), 'methods', {{'prop', 'direct', 'adpf', 'adp'}});
    case 'split'                                      % split predictor (causal LMP average + net-load features), same 36 settings
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        hd = [60 30 7];
        [f, l, d] = ind2sub([6, 2, 3], id);           % 1..36
        cfg = struct('tag', sprintf('split_%s_%d_%d', ft{l}, hd(d), f), 'h', H{f}, 'F', FL(l), ...
            'hist', hd(d), 'methods', {{'split', 'direct'}});
    case 'sens'
        % {label, a, qs, scale, history, eps, refit}
        C = { 'base', 0.40, 0.005, 1.00, 60, 0.32, 30; 'a low', aL(1), 0.005, 1.00, 60, 0.32, 30;
              'a high', aL(2), 0.005, 1.00, 60, 0.32, 30; 'fleet x0.75', 0.40, 0.005, 0.75, 60, 0.32, 30;
              'fleet x1.25', 0.40, 0.005, 1.25, 60, 0.32, 30; 'history 30 d', 0.40, 0.005, 1.00, 30, 0.32, 30;
              'refit 90 d', 0.40, 0.005, 1.00, 60, 0.32, 90; 'eps=0.10', 0.40, 0.005, 1.00, 60, 0.10, 30;
              'eps=0.50', 0.40, 0.005, 1.00, 60, 0.50, 30; 'q+=0.002', 0.40, 0.002, 1.00, 60, 0.32, 30;
              'q+=0.010', 0.40, 0.010, 1.00, 60, 0.32, 30; 'history 7 d', 0.40, 0.005, 1.00, 7, 0.32, 30 };
        [c, l] = ind2sub([12, 2], id);                % 1..24
        cfg = struct('tag', sprintf('sens_%s_%d', ft{l}, c), 'label', C{c,1}, 'a', C{c,2}, 'qs', C{c,3}, ...
            'scale', C{c,4}, 'hist', C{c,5}, 'eps', C{c,6}, 'refit', C{c,7}, 'F', FL(l), ...
            'methods', {{'prop', 'adp', 'myop'}});
    case 'mpcerr'                                     % MPC with forecast errors, F = inf and limited
        [r, l] = ind2sub([4, 2], id); R = [5 10 20 50];  % 1..8
        cfg = struct('tag', sprintf('mpcerr_%s_%d', ft{l}, R(r)), 'F', FL(l), 'methods', {{sprintf('mpc%d', R(r))}});
    case 'mpcfig'                                     % MPC with 10% error for Fig. 2 (1..8)
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        if id <= 6
            cfg = struct('tag', sprintf('mpcfig_fleet_%d', id), 'h', H{id}, 'F', FL(2), 'methods', {{'mpc10'}});
        else
            cfg = struct('tag', sprintf('mpcfig_lim_%g', F2(id-6)), 'F', F2(id-6), 'methods', {{'mpc10'}});
        end
    case 'mpcfit'                                     % MPC with a fitted forecast (1..2)
        cfg = struct('tag', sprintf('mpcfit_%s', ft{id}), 'F', FL(id), 'methods', {{'mpcfit'}});
    case 'trend30'                                    % slow resource on a 30-day mean, line limited (1..3)
        M3 = {{'prop', 'adp'}, {'mpcfit'}, {'mpc20', 'mpc50'}};
        cfg = struct('tag', sprintf('trend30_%d', id), 'F', FL(2), 'trend', 30, 'methods', {M3{id}});
    case 'mpcu'                                       % MPC, forecast = true +- uniform r% of net load (1..12)
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        if id <= 4                                    % Table I: {10, 20}% x {free, lim}
            [r, l] = ind2sub([2, 2], id); R = [10 20];
            cfg = struct('tag', sprintf('mpcu_%s_%d', ft{l}, R(r)), 'F', FL(l), 'methods', {{sprintf('mpcu%d', R(r))}});
        elseif id <= 10                               % Fig. 2 fleets, limited, 10%
            cfg = struct('tag', sprintf('mpcu_fleet_%d', id-4), 'h', H{id-4}, 'F', FL(2), 'methods', {{'mpcu10'}});
        else                                          % Fig. 2 limits 4.2, 2.4 kW, 10%
            cfg = struct('tag', sprintf('mpcu_lim_%g', F2(id-10)), 'F', F2(id-10), 'methods', {{'mpcu10'}});
        end
    case 'mpcfitfig'                                  % fitted-forecast MPC for Fig. 2 (1..8)
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        if id <= 6
            cfg = struct('tag', sprintf('mpcfit_fleet_%d', id), 'h', H{id}, 'F', FL(2), 'methods', {{'mpcfit'}});
        else
            cfg = struct('tag', sprintf('mpcfit_lim_%g', F2(id-6)), 'F', F2(id-6), 'methods', {{'mpcfit'}});
        end
    case 'mpcxgb'                                     % XGBoost-forecast MPC, six fleets x {free, lim} (1..12)
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        [f, l] = ind2sub([6, 2], id);
        cfg = struct('tag', sprintf('mpcxgb_%s_fleet_%d', ft{l}, f), 'h', H{f}, 'F', FL(l), 'methods', {{'mpcxgb'}});
    case 'mpcfitfree'                                 % fitted-forecast MPC, six fleets without limit (1..6)
        H = {[4 30], [1 200], [4 7.8 15.4 30], [1 4 30 200], [4 6 9 13.4 20 30], [1 2.9 8.4 24 69 200]};
        cfg = struct('tag', sprintf('mpcfit_free_fleet_%d', id), 'h', H{id}, 'F', FL(1), 'methods', {{'mpcfit'}});
    case 'eref'                                       % SoC reference 0.5 E, base fleet (1..2)
        cfg = struct('tag', sprintf('eref50_%s', ft{id}), 'F', FL(id), 'eref', 0.5, 'methods', {{'prop', 'adp'}});
    case 'tight'                                      % tighter limits of line 5-6 (1..6)
        Ft = [1.5 1.0 0.5]; [k, m] = ind2sub([3, 2], id); Mt = {{'prop', 'adp', 'blind'}, {'mpcfit'}};
        cfg = struct('tag', sprintf('tight_%g_%d', Ft(k), m), 'F', Ft(k), 'methods', {Mt{m}});
    case 'wide6'                                      % 6-wide fleet under more limits of line 5-6 (1..6)
        Fw = [4.2 2.4 1.5]; [k, m] = ind2sub([3, 2], id); Mw = {{'prop', 'adp'}, {'mpcfit'}};
        cfg = struct('tag', sprintf('wide6_%g_%d', Fw(k), m), 'h', [1 2.9 8.4 24 69 200], 'F', Fw(k), ...
            'methods', {Mw{m}});
    case 'lim'
        cfg = struct('tag', sprintf('lim_%g', F2(id)), 'F', F2(id), 'methods', {{'prop', 'adp', 'mpc', 'myop', 'blind'}});
end
if is14
    cfg.sys = 'ieee14'; cfg.tag = ['i14_' cfg.tag];
    if isfield(cfg, 'a') && cfg.a == 0.40, cfg.a = 0.30; end       % base cost of the units at buses 1, 2
end
if isS
    cfg.tag = ['i14s_' cfg.tag(5:end)]; cfg.trend = 0; cfg.termval = true;
    if strcmp(set, 'main'), cfg.methods = [cfg.methods, {'hindcal', 'mpcfit'}]; cfg.traj = true; end
end
end
