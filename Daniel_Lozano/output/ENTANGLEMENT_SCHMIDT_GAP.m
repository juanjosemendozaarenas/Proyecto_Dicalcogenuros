% Function to read results of the superfermion TNT calculation

function ENTANGLEMENT_SCHMIDT_GAP

clear; clc;

% System and lead parameters
L = 64; % Number of system sites.
J=1;
U=-2;
V=[0,0.5,1,1.5,2,2.5,3,3.5,4];
a=0.1;

Entr=ones(size(V))*-1;%9,1);
Schmidt_Gap=ones(size(V))*-1;%9,1);
disp(size(Schmidt_Gap))
% extra= ['./U' num2str(U) '_L' num2str(L) '/'];
% disp(extra)

for i=1:9
  
%Entang_chi200_GS_FH_NNN_L64_[33_33]_J1_U-2_V1.5_a0.1_chi1000    
fname = ['Entang_chi200_GS_FH_NNN_L' num2str(L) '_[32_32]_J' num2str(J) '_U' num2str(U) '_V' num2str(V(i)) '_a' num2str(a) '_chi1000.mat'];
check = exist(fname);

if check == 2    
    disp([ 'L' num2str(L) '_[32_32]_J' num2str(J) '_U' num2str(U) '_V' num2str(V(i)) '_a' num2str(a) '_chi1000.mat']);
    load(fname); % Load file
    Entr(i)=Entropy(end);
    Schmidt_Gap(i)=Schmidt_coeff(1).^2-Schmidt_coeff(2).^2;
    disp( ["Entropy=" num2str(Entropy(end)) ] );
    disp(["Schmidt gap= " num2str(Schmidt_Gap(i))]);
end

end

hold on 
P1=plot(V,Entr,"-o");
L1="Entropy";
P2=plot(V,Schmidt_Gap,"-o");
L2="Schmidt gap";
grid on
legend([P1,P2],[L1,L2]);