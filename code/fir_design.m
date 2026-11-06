function out = fir_design(ac,q,c,kap,dt,L,twoSided,rho)
% Exact fixed-price, zero-DC FIR synthesis for a stated tap support.
% J includes rho*||theta||^2. Q is its quadratic coefficient, not Hessian.
persistent cachedAc cachedDt cachedM cachedR cachedG
if twoSided, m=2*L-1; off=-(L-1):(L-1); else, m=L; off=0:L-1; end
ac=ac(:); N=numel(c);
assert(numel(ac)>=m,'Insufficient autocovariance lags.');
if isempty(cachedM) || cachedM~=m || cachedDt~=dt || ~isequal(cachedAc,ac)
    R=toeplitz(ac(1:m)); C=tril(ones(m)); C=C(1:end-1,:);
    G=dt^2*(C'*toeplitz(ac(1:m-1))*C);
    cachedAc=ac; cachedDt=dt; cachedM=m; cachedR=R; cachedG=G;
else
    R=cachedR; G=cachedG;
end
Q=kron(q*ones(N)+diag(c),R)+kron(diag(kap),G)+rho*eye(N*m);
Q=(Q+Q')/2;
g=repmat(q*ac(abs(off)+1),N,1);
F=kron(eye(N),ones(1,m));
D=decomposition(Q,'chol'); U=D\[g,F'];
lambda=(F*U(:,2:end))\(F*U(:,1));
theta=U(:,1)-U(:,2:end)*lambda;
out.theta=reshape(theta,m,N);
out.J=q*ac(1)-2*g'*theta+theta'*Q*theta;
out.sigP=sqrt(max(sum(out.theta.*(R*out.theta),1),0));
out.sigE=sqrt(max(sum(out.theta.*(G*out.theta),1),0));
out.Q=Q;
out.kktResidual=max(norm(Q*theta+F'*lambda-g,inf),norm(F*theta,inf));
end
