function [sx,sy,sz,sp,sm] = tntMatPauliOps(s,N)
%% tntMatPauliOps Creates the Pauli matrices for N-species spin operators with a spin given by TwoS

% z component of spin
m = s:-1:-s;

% Build sparse matrix version of basic spin operators :
sz = 2*diag(m);
sp = diag(sqrt((s-m(2:end)).*(s+m(2:end)+1)),1);
sm = diag(sqrt((s+m(1:(end-1))).*(s-m(1:(end-1))+1)),-1);

sx = sp+sm;
sy = -1i*(sp-sm);

% Construct spin operators for each spin in the full Hilbert space :
sx = tntMatExpandBasis(sx,N);
sy = tntMatExpandBasis(sy,N);
sz = tntMatExpandBasis(sz,N);
sp = tntMatExpandBasis(sp,N);
sm = tntMatExpandBasis(sm,N);
