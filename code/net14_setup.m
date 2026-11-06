function net = net14_setup(F)
% NET14_SETUP  IEEE 14-bus system (MATPOWER case14) scaled to the CoSSMic
% case. DC power flow from the branch reactances, slack bus 1. The eleven
% load buses carry net-load deviations proportional to their case14 loads
% (net14_loads); the five case14 generators keep their buses, with the
% case14 costs scaled so that the units at buses 1 and 2 cost 0.30 $/kWh
% and those at buses 3, 6, 8 cost 0.60 $/kWh; quadratic terms keep the
% case14 ratios with b = 0.005 $/(kW^2 h) at bus 1. Load buses can curtail.
%   F  line limits (20 x 1, kW); default: no limits.
% Data transcribed from MATPOWER data/case14.m (IEEE CDF, UW archive).
br = [ 1 2 0.05917;  1 5 0.22304;  2 3 0.19797;  2 4 0.17632;  2 5 0.17388;
       3 4 0.17103;  4 5 0.04211;  4 7 0.20912;  4 9 0.55618;  5 6 0.25202;
       6 11 0.1989;  6 12 0.25581; 6 13 0.13027; 7 8 0.17615;  7 9 0.11001;
       9 10 0.0845;  9 14 0.27038; 10 11 0.19207; 12 13 0.19988; 13 14 0.34802 ];
Pd = [0 21.7 94.2 47.8 7.6 11.2 0 0 29.5 9 3.5 6.1 13.5 14.9];      % MW, case14
c1 = [20 20 40 40 40]; c2 = [0.0430292599 0.25 0.01 0.01 0.01];     % case14 gencost
nb = 14; L = size(br, 1);
Bf = zeros(L, nb); Bbus = zeros(nb);
for l = 1:L
    f = br(l,1); k = br(l,2); b = 1/br(l,3);
    Bf(l,f) = b; Bf(l,k) = -b;
    Bbus(f,f) = Bbus(f,f) + b; Bbus(k,k) = Bbus(k,k) + b;
    Bbus(f,k) = Bbus(f,k) - b; Bbus(k,f) = Bbus(k,f) - b;
end
PTDF = zeros(L, nb);
PTDF(:, 2:nb) = Bf(:, 2:nb) / Bbus(2:nb, 2:nb);    % injection at n, withdrawal at slack bus 1
net.lines = br; net.nb = nb; net.L = L; net.PTDF = PTDF;
net.names = arrayfun(@(n) sprintf('%d', n), 1:nb, 'UniformOutput', false);
net.loadBus = find(Pd > 0);                         % 2 3 4 5 6 9 10 11 12 13 14
net.Pd = Pd(net.loadBus);
net.unitBus = [3 4 9 14];                           % 1 h, 4 h, 30 h, 200 h (netrun may override)
net.genBus  = [1 2 3 6 8];
net.genA    = 0.015 * c1;                           % 0.30, 0.30, 0.60, 0.60, 0.60 $/kWh
net.genQ    = 0.005 * c2 / c2(1);                   % $/(kW^2 h)
net.nGen    = 2;                                    % units at buses 1, 2 are the main generators
net.curtBus = net.loadBus;
net.qc      = 0.002;
if nargin < 1 || isempty(F), F = inf(L, 1); end
net.F = F(:);
end
