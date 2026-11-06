function F = price_features(x, hrs, memH)
% PRICE_FEATURES  Causal regressors at each step t (after x_t is observed):
% hourly means of the last memH hours of net load (newest first; zero
% where the history is shorter) and the first two daily harmonics of the
% hour of day hrs(t).
n = numel(x); c = [0; cumsum(x(:))];
F = zeros(n, memH + 4);
for hh = 1:memH
    b = (1:n)' - (hh-1)*4; a = b - 3;
    ok = a >= 1;
    F(ok, hh) = (c(b(ok)+1) - c(a(ok)))/4;
end
h = hrs(:);
F(:, memH+1:memH+4) = [sin(2*pi*h/24), cos(2*pi*h/24), sin(4*pi*h/24), cos(4*pi*h/24)];
end
