function out = refit_values(xs, hrs, par, nTrain0, blkDays, variants)
% REFIT_VALUES  Monthly refit on all past data (expanding window). At the
% start of each block of the test period, one hindsight QP on xs(1:t0) is
% shared by all price-form variants and by the ADP value fit.
%   variants{v}.pol    prices (kappa, ceff) of variant v
%   variants{v}.opts   options for design_price
%   variants{v}.recal  (optional, v26) if true, recalibrate kappa and ceff
%                      with synth_policy on xs(1:t0) at every refit
%   variants{v}.calMin (optional, v26b) days of data the calibration needs.
%                      Before that, the variant uses the duration prior
%                      kappa_i = c_i/h_i^2 (tau_i = h_i), the starting point
%                      of synth_policy; it calibrates once at the first
%                      refit with at least calMin days and keeps the result.
%   variants{v}.minHist (optional, v26c) days of data the predictor needs
%                      (its lags span 8 days); blocks that start with less
%                      use constant values (window mean of the SoC duals).
% Returns out.V{v} (T x N unit values), out.fits{b} (ADP values), and
% out.pol{v,b} (prices used by variant v in block b).
spd = 96; T = numel(xs) - nTrain0; blk = blkDays*spd; nb = ceil(T/blk);
nv = numel(variants);
out.V = cell(1, nv);
for v = 1:nv, out.V{v} = zeros(T, par.N); end
out.fits = cell(1, nb); out.blk = blk; out.cal = cell(nv, nb); out.pol = cell(nv, nb);
prior = struct('kappa', max(par.c(:)'./par.h(:)'.^2, 1e-9), 'ceff', par.c(:)');
polKeep = cell(1, nv);
for b = 1:nb
    b0 = (b-1)*blk + 1; b1 = min(T, b*blk); t0 = nTrain0 + b0 - 1;
    w = xs(1:t0);
    [~, Ht] = evalc('method_hindsight(w, par)');
    out.fits{b} = fit_values(w, par, Ht);
    polw = [];
    for v = 1:nv
        pol = variants{v}.pol;
        if isfield(variants{v}, 'recal') && variants{v}.recal
            if isempty(polw), [~, polw] = evalc('synth_policy(w, par)'); end
            pol = polw;
        elseif isfield(variants{v}, 'calMin')
            if isempty(polKeep{v})
                if t0 >= variants{v}.calMin*spd
                    if isempty(polw), [~, polw] = evalc('synth_policy(w, par)'); end
                    polKeep{v} = polw; pol = polw;
                else
                    pol = prior;
                end
            else
                pol = polKeep{v};
            end
        end
        o = variants{v}.opts; o.Ht = Ht;
        if isfield(variants{v}, 'minHist') && t0 < variants{v}.minHist*spd, o.constant = true; end
        pp = design_price(w, hrs(1:t0), par, pol, o);
        out.V{v}(b0:b1, :) = predict_values(xs(t0+1:nTrain0+b1), w, hrs(1:nTrain0+b1), pp);
        out.cal{v, b} = pp.cal; out.pol{v, b} = pol;
    end
end
end
