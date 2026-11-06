function V = predict_values(x, xhist, hrs, pp)
% PREDICT_VALUES  Unit values of the price-form policy for the samples x,
% given the preceding history xhist and the hours hrs of [xhist; x].
% With pp.multi (network), x and xhist hold one column per bus.
if isfield(pp, 'multi') && pp.multi
    xa = [xhist; x];
    F = price_features_multi(xa, hrs, pp.memH);
    F = F(size(xhist,1)+1:end, :);
    V = [ones(size(x,1),1), (F - pp.mu)./pp.sd] * pp.B;
    return
end
xa = [xhist(:); x(:)];
F = price_features(xa, hrs, pp.memH);
F = F(numel(xhist)+1:end, :);
if isfield(pp, 'useHarm') && ~pp.useHarm, F = F(:, 1:pp.memH); end
V = [ones(numel(x),1), (F - pp.mu)./pp.sd] * pp.B;
end
