function F = price_features_multi(X, hrs, memH)
% PRICE_FEATURES_MULTI  price_features on several causal series (one per
% bus): the hourly lags of each column, then the four daily harmonics once.
k = size(X, 2); n = size(X, 1);
F = zeros(n, k*memH + 4);
for c = 1:k
    Fc = price_features(X(:,c), hrs, memH);
    F(:, (c-1)*memH + (1:memH)) = Fc(:, 1:memH);
    if c == 1, F(:, k*memH + (1:4)) = Fc(:, memH+1:memH+4); end
end
end
