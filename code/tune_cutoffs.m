function [wBest, tab] = tune_cutoffs(xval, par, pol0)
% TUNE_CUTOFFS  Select the three moving-average windows of the fixed-cutoff
% baseline on the validation segment (same budget as the proposed gain):
% first cutoff 0.5-16 h, second 6-96 h, third 2-28 days, with step 3.
W1 = [2 4 8 16 32 64]; W2 = [24 48 96 192 384]; W3 = [192 288 672 1344 2688];
best = inf; wBest = [4 48 288]; tab = [];
for a = W1
    for b = W2(W2 > a)
        for c = W3(W3 > b)
            ph = pol0; ph.maWin = [a b c]; ph.completion = true;
            ev = evalc('R = method_filter_heuristic(xval, par, ph);'); %#ok<NASGU>
            tab(end+1,:) = [a b c R.cost]; %#ok<AGROW>
            if R.cost < best, best = R.cost; wBest = [a b c]; end
        end
    end
end
fprintf('[cutoffs] selected [%s] steps (%.1f $ on validation)\n', ...
    sprintf('%d ', wBest), best);
end
