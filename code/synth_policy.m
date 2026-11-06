function pol = synth_policy(xtrain, par)
% Fixed-class causal FIR synthesis; statistics use the supplied past only.
% Variance targets are design proxies, not guarantees for the clipped policy.
L=768; if isfield(par,'firLength'), L=par.firLength; end
rho=5e-10;
if isfield(par,'design_z')
 z=par.design_z;
else
 % Per-limit design target: each of the power and energy limits separately
 % (joint coverage follows only from a union bound, 1-2*eps).
 z=sqrt(2)*erfinv(1-par.eps_risk);
end
mu=mean(xtrain); xz=xtrain(:)-mu; n=numel(xz);
assert(n>=2 && L>=2,'At least two samples and taps are required.');
nfft=2^nextpow2(max(2*n,4*L));
ac=real(ifft(abs(fft(xz,nfft)).^2))/n;
ac=ac(1:2*L); ac(n+1:end)=0;
kap=max(par.c./par.h.^2,1e-9); ceff=par.c;
head=min(par.eref,par.E-par.eref);
assert(all(head>0),'Reference SoC must be inside the energy limits.');
converged=false;
for it=1:80
 d=fir_design(ac,par.q,ceff,kap,par.dt,L,false,rho);
 ve=z*d.sigE./head; vp=z*d.sigP./par.P;
 tightE=all(ve>=0.85 | kap<=1.001e-9);
 tightP=all(vp>=0.85 | ceff<=par.c*1.001);
 if max([ve,vp])<=1.001 && tightE && tightP
  converged=true; break;
 end
 kap=max(kap.*max(min(ve.^2,25),0.04),1e-9);
 ceff=max(ceff.*max(min(vp.^2,25),0.20),par.c);
end
% Resolve at the final stored prices, including after the iteration cap.
for retry=1:25
 d=fir_design(ac,par.q,ceff,kap,par.dt,L,false,rho);
 ve=z*d.sigE./head; vp=z*d.sigP./par.P;
 if max([ve,vp])<=1.005, converged=true; break; end
 % Uniform price growth is a safe fallback toward zero dispatch.
 scale=max(1.2,max([ve,vp])^2);
 kap=kap*scale; ceff=ceff*scale;
end
assert(converged,'Variance calibration failed: max ratio %.4f',max([ve,vp]));
pol=struct('theta',d.theta,'L',L,'kappa',kap,'ceff',ceff,'mu',mu,...
 'sigE',d.sigE,'sigP',d.sigP,'z',z,'rho',rho,'Jc_pred',d.J,...
 'gamma',par.gamma,'calibration_converged',converged,...
 'calibration_iterations',it,'variance_ratios',[ve;vp],...
 'history',xtrain(max(1,end-L+1):end));
clear d
b=fir_design(ac,par.q,ceff,kap,par.dt,L,true,rho);
pol.thetaTwo=b.theta; pol.Jnc_pred=b.J;
pol.lqGap=pol.Jc_pred-pol.Jnc_pred;
pol.dGap_pred=pol.lqGap; % compatibility field: now the same regularized LQ objective
embedded=[zeros(L-1,par.N);pol.theta]; diff=embedded(:)-b.theta(:);
pol.gap_identity_error=abs(pol.lqGap-diff'*b.Q*diff);
assert(pol.lqGap>=-1e-7 && pol.gap_identity_error<1e-5*max(1,abs(pol.Jc_pred)),...
 'Causal/two-sided FIR gap identity failed.');
pol.ac=ac;
pol.f=logspace(-4,log10(2),400); w=2*pi*pol.f;
ai=ceff(:)+kap(:)./w.^2; gi=1./ai;
pol.H=gi./(1/par.q+sum(gi,1));
pol.K=(exp(-1j*(w(:)*par.dt)*(0:L-1))*pol.theta).';
fprintf('[synth] L=%d, iterations=%d, target ratio %.4f, LQ gap %.6f per hour, identity error %.2g\n',...
 L,it,max([ve,vp]),pol.lqGap,pol.gap_identity_error);
end
