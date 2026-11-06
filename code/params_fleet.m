function par = params_fleet(x, dt)
% PARAMS_FLEET  Four-duration storage fleet scaled to the net-load statistics.
% Units sorted short -> long duration; c increasing, kappa (calibrated later)
% decreasing. Costs in $/(kW^2 h); q chosen so shedding marginal cost at
% s = 10 kW matches a VoLL-like 10 $/kWh.

sx = std(x);

par.dt = dt;
par.N  = 4;
par.name = {'BESS-P (1 h)', 'BESS-E (4 h)', 'Flow (30 h)', 'H2 (200 h)'};

% power ratings relative to net-load std; energy = P * duration.
% Deliberately scarce: the fleet cannot absorb everything at once, so both
% power and energy limits bind and intertemporal allocation matters
% (with abundant capacity, myopic merit-order dispatch is already optimal).
f       = [0.8, 0.6, 0.45, 0.35];
h       = [1, 4, 30, 200];               % durations [h]
par.P   = round(f * sx, 1);              % kW
par.E   = par.P .* h;                    % kWh
par.h   = h;

par.etac = [0.96, 0.95, 0.85, 0.65];
par.etad = [0.96, 0.95, 0.85, 0.55];

% quadratic cycling-cost coefficients [$ / (kW^2 h)] (convex aging surrogate)
par.c0  = [0.001, 0.003, 0.010, 0.030];
% first-order loss folding: c_i <- c_i + ctilde * (1 - etac*etad)
par.c   = par.c0 + 0.01 * (1 - par.etac .* par.etad);

% residual imbalance economics: s>0 is served by the backup genset at a
% linear fuel cost (plus small convex wear curvature); s<0 is curtailed
% surplus PV at near-zero cost. Storage value = displaced fuel - cycling.
par.a_fuel = 0.40;                       % $/kWh genset marginal fuel cost
% small off-grid genset: efficiency degrades sharply with loading, so the
% fuel cost is strongly convex (marginal cost 4x the base at s = +30 kW).
% This is what makes intertemporal allocation valuable: peaks must be
% pre-positioned for, and misallocated bands are genuinely punished.
par.qs  = 0.020;                         % $/(kW^2 h) genset convexity
par.qc  = 0.002;                         % $/(kW^2 h) curtailment curvature
% symmetric quadratic surrogate for synthesis: matches the marginal fuel
% cost at one standard deviation, 2*q*sx = a + 2*qs*sx
par.q   = par.a_fuel / (2*sx) + par.qs;
% chance level for calibration: deliberately loose (z ~ 1). Violations are
% benign -- the projection layer clips them -- while a tight level makes the
% policy too timid to use its capacity (the LQ policy rides E/2 +- z*sigma,
% whereas the constrained optimum parks at the SoC bounds).
par.eps_risk = 0.32;

% SoC anti-windup gains: below each unit's band edge (1/h_i per hour)
par.gamma = 0.1 ./ h;                    % kW per kWh of SoC deviation

% asymmetric imbalance costs shift the optimal SoC set-point above midpoint
% (shedding dear, curtailment cheap): bias stored energy upward.
par.eref = 0.6 * par.E;
par.e0   = par.eref;                     % start at reference SoC
end
