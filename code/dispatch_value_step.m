function [p,s,e] = dispatch_value_step(x,e_prev,value,par)
% Exact one-step dispatch with a nonnegative linear stored-energy value.
% Charge/discharge bounds are inside the common-price balance equation.
pmax=min(par.P,par.etad.*e_prev/par.dt);
pmin=-min(par.P,(par.E-e_prev)./(par.etac*par.dt));
value=max(value,0);
hi=max([1,par.a_fuel+2*par.qs*abs(x),value./par.etad+2*par.c0.*par.P]);
lo=-max([1,2*par.qc*abs(x),2*par.c0.*par.P]);
for it=1:52
    mu=(lo+hi)/2;
    discharge=max((mu-value./par.etad)./(2*par.c0),0);
    charge=min((mu-value.*par.etac)./(2*par.c0),0);
    req=min(pmax,max(pmin,discharge+charge));
    if mu>par.a_fuel
        sr=(mu-par.a_fuel)/(2*par.qs);
    elseif mu<0
        sr=mu/(2*par.qc);
    else
        sr=0;
    end
    if sum(req)+sr<x, lo=mu; else, hi=mu; end
end
[p,s,e]=util_clip(req,x,e_prev,par);
end
