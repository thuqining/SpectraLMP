function [p, s, e] = util_clip(ptil, xt, e_prev, par)
% UTIL_CLIP  Physical feasibility only: each unit's request is clipped
% independently to its power/SoC limits; the genset/curtailment absorbs
% the residual. No cross-unit reallocation -- coordination is each
% policy's own responsibility (the proposed algorithm performs its
% merit-order reallocation as its own step 3, see util_project).

dt = par.dt;
pmax = min(par.P,  par.etad .* e_prev / dt);
pmin = -min(par.P, (par.E - e_prev) ./ (par.etac * dt));
p = min(pmax, max(pmin, ptil));
s = xt - sum(p);

e = e_prev;
dis = p > 0;
e(dis)  = e(dis)  - p(dis) ./ par.etad(dis) * dt;
e(~dis) = e(~dis) - p(~dis) .* par.etac(~dis) * dt;
e = min(max(e, 0), par.E);
end
