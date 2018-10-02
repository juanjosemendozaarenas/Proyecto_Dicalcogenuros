function [HAMILTONIAN,NUMERO] = FERMI_HUBBARD_FUNCTION(C_d_do,C_do,N_do,C_d_up,C_up,N_up,P,A_id,J,U,V,dim,N,PBC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ham=zeros(dim,dim);
Numero=zeros(dim,dim);

CdC=kron(C_d_do*P,C_do)+kron(C_d_up*P,C_up);%primer termino de hopping

CCd=kron(P*C_do,C_d_do)+kron(P*C_up,C_d_up);%segundo termino de hopping

N_T=N_do+N_up;%Numero de fermiones por sitio

N_T_2=kron(N_T,N_T);%interaccion densidad densidad

for i=(1:N-1)
    ID1=A_id{i};
    ID2=A_id{N-i};
    Ham=Ham-J*kron(ID1,kron(CdC+CCd,ID2))+V*kron(ID1,kron(N_T_2,ID2));
end

N_up_do=N_up*N_do;
%N_do_up=N_do*N_up;

for i=(1:N)
   ID1=A_id{i};
   ID2=A_id{N-i+1};
   Ham=Ham+U*(kron(ID1,kron(N_up_do,ID2))) ;
   Numero=Numero+kron(ID1,kron(N_T,ID2));
end

P_in=P;
for i=(1:N-3)
P_in=kron(P_in,P);
end

Ham=Ham+PBC*(-J*(kron(P*(C_do),kron(P_in,C_d_do)) ...
                +kron(P*(C_up),kron(P_in,C_d_up))...
                +kron((C_d_do)*P,kron(P_in,C_do)) ...
                +kron((C_d_up)*P,kron(P_in,C_up))) ...
            +V*(kron((N_T),kron(A_id{N-1},N_T) ) )  );
HAMILTONIAN=Ham;
NUMERO=Numero;



end

