function net = net5_setup(F)
% NET5_SETUP  PJM 5-bus test network (buses A-E = 1-5) scaled to the
% CoSSMic case. DC power flow with slack bus A. Net loads of the three
% buildings sit at the load buses B, C, D; storage units (1 h, 4 h, 30 h,
% 200 h) sit at B, C, D, E; a cheap generator sits at E and the
% case-study generator at A; curtailment is possible at each load bus.
%   F  line limits (6 x 1, kW); default: no limits.
lines = [1 2 0.0281; 1 4 0.0304; 1 5 0.0064; 2 3 0.0108; 3 4 0.0297; 4 5 0.0297];
nb = 5; L = size(lines, 1);
Bf = zeros(L, nb); Bbus = zeros(nb);
for l = 1:L
    f = lines(l,1); k = lines(l,2); b = 1/lines(l,3);
    Bf(l,f) = b; Bf(l,k) = -b;
    Bbus(f,f) = Bbus(f,f) + b; Bbus(k,k) = Bbus(k,k) + b;
    Bbus(f,k) = Bbus(f,k) - b; Bbus(k,f) = Bbus(k,f) - b;
end
PTDF = zeros(L, nb);
PTDF(:, 2:nb) = Bf(:, 2:nb) / Bbus(2:nb, 2:nb);    % injection at n, withdrawal at slack A
net.lines = lines; net.nb = nb; net.L = L; net.PTDF = PTDF;
net.names = {'A','B','C','D','E'};
net.loadBus = [2 3 4];                  % ind1, res5, res6
net.unitBus = [2 3 4 5];                % 1 h, 4 h, 30 h, 200 h
net.genBus  = [1 5 2 3 4];              % main generators at A and E; local peakers at B, C, D
net.genA    = [0.40 0.30 0.60 0.60 0.60];  % $/kWh
net.genQ    = [0.005 0.005 0.005 0.005 0.005];  % $/(kW^2 h)
net.nGen    = 2;                        % the first nGen entries are the main generators
net.curtBus = [2 3 4];                  % curtailment at the load buses
net.qc      = 0.002;                    % $/(kW^2 h)
if nargin < 1 || isempty(F), F = inf(L, 1); end
net.F = F(:);
end
