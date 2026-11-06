function c = net5_cost(pd, pc, sg, sc, par, net)
% NET5_COST  Operating cost (1a) of a network dispatch, $.
dt = par.dt;
c = dt * ( sum(par.c0(:) .* (pd + pc).^2, 'all') ...
         + sum(net.genA(:) .* sg + net.genQ(:) .* sg.^2, 'all') ...
         + net.qc * sum(sc.^2, 'all') );
% Optional terminal value: energy drawn from (added to) the initial stock is
% charged (credited) at the lowest producing generator cost times etad, the
% value Step 2 of Algorithm 1 gives to the final energy.
if isfield(par, 'termval') && par.termval
    eT = par.e0(:) + dt * sum(par.etac(:) .* pc - pd ./ par.etad(:), 2);
    c = c + par.a_fuel * sum(par.etad(:) .* (par.e0(:) - eT));
end
end
