function [n,c,cd,Pc,cdP,sz,sz2,dbl] = tntMatFermionOps
%% tntMatFermionOps Creates operators of spin-1/2 fermions (single species)

nup = diag([0,1,0,1]);
ndn = diag([0,0,1,1]);
cup = [0,1,0,0;
       0,0,0,0;
       0,0,0,1;
       0,0,0,0];
cdup = cup';
cdn = [0,0,1,0;
       0,0,0,-1;
       0,0,0,0;
       0,0,0,0];
cddn = cdn';
sz = nup - ndn;
sz2 = sz*sz;
dbl = nup*ndn;
P = (-1)^(nup+ndn); 

% Turn into cells, to obtain correct TNT form:

% Separated operators of both species
n = {ndn;nup};
c = {cdn;cup};
cd = {cddn;cdup};

% Mixed operators of both species
sz = tntMatExpandBasis(sz,1);
sz2 = tntMatExpandBasis(sz2,1);
dbl = tntMatExpandBasis(dbl,1);

% Operators for kinetic term
Pc = {P*cdn;P*cup}; % Pcdn and Pcup
cdP = {cddn*P;cdup*P}; % cddnP and cdupP