% Function to read results of the superfermion TNT calculation

function GAP_C_S

clear; clc;

% System and lead parameters
L = 64; % Number of system sites.
J=1;
U=-2;
V=[0,0.5,1,1.5,2,2.5,3,3.5,4];
DS=ones(size(V))*-1;
DC=ones(size(V))*-1;
a=0.1;
chi1=1000;
chi2=2000;

Number1=[31,32,33,33,31];
Number2=[31,32,33,31,33];

% Search for the V files which have 4 output files
for j=1:9
file=['V' num2str(V(j))];
disp(file);
path(path,file); % Add path for the V file

Energy=zeros(5,1);
counter=0;

% Searching for output files
for i=1:5
    
fname_1 = ['V' num2str(V(j)) '/' 'GS_FH_NNN_L' num2str(L) '_[' num2str(Number1(i)) '_' num2str(Number2(i)) ']_J' num2str(J) '_U' num2str(U) '_V' num2str(V(j)) '_a' num2str(a) '_chi' num2str(chi1) '.mat'];
check1 = exist(fname_1);
fname_2 = ['V' num2str(V(j)) '/' 'GS_FH_NNN_L' num2str(L) '_[' num2str(Number1(i)) '_' num2str(Number2(i)) ']_J' num2str(J) '_U' num2str(U) '_V' num2str(V(j)) '_a' num2str(a) '_chi' num2str(chi2) '.mat'];
check2 = exist(fname_2);

if check1 == 2    
    %disp([ 'L' num2str(L) '_[' num2str(Number1(i)) '_' num2str(Number2(i)) ']_J' num2str(J) '_U' num2str(U) '_V' num2str(V(j)) '_a' num2str(a)]);
    counter=counter+1;
    load(fname_1); % Load file
    Energy(i)=E(end);
    disp(E(end))
end

if check2 == 2    
    %disp([ 'L' num2str(L) '_[' num2str(Number1(i)) '_' num2str(Number2(i)) ']_J' num2str(J) '_U' num2str(U) '_V' num2str(V(j)) '_a' num2str(a)]);
    counter=counter+1;
    load(fname_2); % Load file
    Energy(i)=E(end);
    disp(E(end))
end
end

% If the output files are enough, the gaps are calculated
    if counter==4
        DC(j)= 0.5 * ( Energy(1)+ Energy(3) -2*Energy(2) );
        
% In case one doesn't work    
    if Energy(4)==0
             DS(j)=Energy(5)-Energy(2);
    else
             DS(j)=Energy(4)-Energy(2);
    end

    
    end


end
hold on
P1=plot(V,DS,"-o"); M1="Spin Gap";
P2=plot(V,DC,"-o"); M2="Charge Gap";
legend([P1,P2],[M1,M2]);
