clc
clear all
close all
DSMCout=importdata('2d_out.tec');



data(:,1)=DSMCout(:,1);
data(:,2)=DSMCout(:,2);
data(:,3)=DSMCout(:,3);
data(:,4)=DSMCout(:,4);
data(:,5)=DSMCout(:,7);
data(:,6)=DSMCout(:,5).*DSMCout(:,15);
data(:,7)=DSMCout(:,5).*DSMCout(:,21);
data(:,8)=DSMCout(:,5).*DSMCout(:,27);
data(:,9)=DSMCout(:,5).*DSMCout(:,33);
data(:,10)=DSMCout(:,5).*DSMCout(:,39);
data(:,11)=DSMCout(:,10);

%% constants and properties of gases
k=1.38064852e-23; %Boltzmann constant
R=288.26; % gas constant for air
diam_air=4.467; % in Ansgtrom
diam_N2=4.467;  % in Ansgtrom 
diam_NO=4.334;  % in Ansgtrom
Molar_air=29;  % molar mass in g/mol 
Molar_N2=28;   % molar mass in g/mol 
Molar_NO=30;   % molar mass in g/mol 
rho_lim_N2A=1e-12; % lower limit for the calculation of Diff coeffs
rho_lim_NOA=1e-15; % lower limit for the calculation of Diff coeffs
R_N=593.89;
R_O=519.65;
R_N2=296.95;
R_O2=259.83;
R_NO=277.15;
Avocado=6.022e+23;

%% geometry and grid variables
Radius=0.015875; % radius of the cylinder
nCellx=284; %number of points in x direction
nCelly=250; %number of points in y direction
nTime=1e+5; %total number of time steps
dt=1e-9; %time step
dx=1e-4; %length betweeen points in x direction
dy=1e-4; %length betweeen points in y direction

%% Calculate rate coefficients, diffusion coefficients and initialization of matrices

% Initialize densities of N2A and NOA
rho_N=zeros(nCellx+1,nCelly+2);
rho_O=zeros(nCellx+1,nCelly+2);
rho_N2=zeros(nCellx+1,nCelly+2);
rho_O2=zeros(nCellx+1,nCelly+2);
rho_NO=zeros(nCellx+1,nCelly+2);
rho_NOA_old=zeros(nCellx+1,nCelly+2);
rho_N2A_old=zeros(nCellx+1,nCelly+2);
rho_NOA_new=zeros(nCellx+1,nCelly+2);
rho_N2A_new=zeros(nCellx+1,nCelly+2);
N_N=zeros(nCellx+1,nCelly+2);
N_O=zeros(nCellx+1,nCelly+2);
N_N2=zeros(nCellx+1,nCelly+2);
N_O2=zeros(nCellx+1,nCelly+2);
N_NO=zeros(nCellx+1,nCelly+2);
N_N2A_old=zeros(nCellx+1,nCelly+2);
N_NOA_old=zeros(nCellx+1,nCelly+2);
N_N2A_new=zeros(nCellx+1,nCelly+2);
N_NOA_new=zeros(nCellx+1,nCelly+2);
D_NOA=zeros(nCellx+1,nCelly+2);
D_N2A=zeros(nCellx+1,nCelly+2);
S_NOA=zeros(nCellx+1,nCelly+2);
S_N2A=zeros(nCellx+1,nCelly+2);
conv_NOA=zeros(nCellx+1,nCelly+2);
conv_N2A=zeros(nCellx+1,nCelly+2);
diff_NOA=zeros(nCellx+1,nCelly+2);
diff_N2A=zeros(nCellx+1,nCelly+2);
NOA70=zeros(nTime,1);
N2A70=zeros(nTime,1);
SNOA70=zeros(nTime,1);
SN2A70=zeros(nTime,1);
convNOA70=zeros(nTime,1);
convN2A70=zeros(nTime,1);
diffNOA70=zeros(nTime,1);
diffN2A70=zeros(nTime,1);
NOA=zeros(nTime,1);
N2A=zeros(nTime,1);




% major quantities
u=zeros(nCellx+1,nCelly+2);
v=zeros(nCellx+1,nCelly+2);
coordX=zeros(nCellx+1,nCelly+2);
coordY=zeros(nCellx+1,nCelly+2);




% Get Number densities of major components and temperature and pressure
% field
P=data(:,11);
Patm=P*1e-5;
T=data(:,5);
NN_N=data(:,6);
NN_O=data(:,7);
NN_N2=data(:,8);
NN_O2=data(:,9);
NN_NO=data(:,10);

% convert to densities in kg/m^3
% 
r_N=NN_N*k/R_N;
r_O=NN_O*k/R_O;
r_N2=NN_N2*k/R_N2;
r_O2=NN_O2*k/R_O2;
r_NO=NN_NO*k/R_NO;



%forward rate coefficients for eqns on Table 3 respectively

Aaf=[8.64e-23;2.48e-43;8.14e-23;1.27e-39;6.64e-17;...
     1.53e-8 ;1.83e-18;4.88e-17;2.32e-18;2.91e-17;...
     1.15e-17;9.24e-18;3.06e-17;1.69e-17;8.64e-18;...
     4.7e+6;7.7e-2];
nf=[-0.35;-1.24;-0.35;-1.60;0.00;-2.23;-0.50;0.00;0.00;0.00;...
     0.50;0.50;0.50;0.50;0.50;0.00;0.00];
 
Eaf=[0.00;0.00;0.00;0.00;0.00;1e-18;9.88e-19;0.00;0.00;0.00;...
     9.70e-20;0.00;0.00;0.00;0.00;0.00;0.00];

Aab=[2.09e+8;5.86e-13;2.33e+9;3.64e-8;1.84e-17;1.31e-10;1.56e-20;...
     5.72e-15;2.73e-16;2.87e-15;4.98e-16;2.87e-15;1.30e-14;...
     5.60e-15;3.64e-15];
 
 
nb=[-0.265;-1.152;-0.692;-1.942;0.471;-1.942;-0.212;-0.289;-0.289;...
    -0.270;0.080;-0.2118;-0.2587;-0.2288;-0.2587];
 
Eab=[1.47e-19;1.47e-19;5.70e-19;5.70e-19;1.01e-19;-1.90e-20;...
     -3.10e-20;1.02e-18;1.02e-18;1.02e-18;1.04e-18;9.21e-19;...
     9.18e-19;9.22e-19;9.18e-19];



%calculate forward and backward coefficients

% initialize rate coefficients
kk1f=zeros(nCellx*nCelly,1);
kk2f=zeros(nCellx*nCelly,1);
kk3f=zeros(nCellx*nCelly,1);
kk4f=zeros(nCellx*nCelly,1);
kk5f=zeros(nCellx*nCelly,1);
kk6f=zeros(nCellx*nCelly,1);
kk7f=zeros(nCellx*nCelly,1);
kk8f=zeros(nCellx*nCelly,1);
kk9f=zeros(nCellx*nCelly,1);
kk10f=zeros(nCellx*nCelly,1);
kk11f=zeros(nCellx*nCelly,1);  
kk12f=zeros(nCellx*nCelly,1);
kk13f=zeros(nCellx*nCelly,1);
kk14f=zeros(nCellx*nCelly,1);
kk15f=zeros(nCellx*nCelly,1);
kk16f=zeros(nCellx*nCelly,1);   
kk17f=zeros(nCellx*nCelly,1);   

kk1b=zeros(nCellx*nCelly,1);
kk2b=zeros(nCellx*nCelly,1);
kk3b=zeros(nCellx*nCelly,1);
kk4b=zeros(nCellx*nCelly,1);
kk5b=zeros(nCellx*nCelly,1);  
kk6b=zeros(nCellx*nCelly,1);
kk7b=zeros(nCellx*nCelly,1);
kk8b=zeros(nCellx*nCelly,1);
kk9b=zeros(nCellx*nCelly,1);
kk10b=zeros(nCellx*nCelly,1);  
kk11b=zeros(nCellx*nCelly,1); 
kk12b=zeros(nCellx*nCelly,1);
kk13b=zeros(nCellx*nCelly,1);
kk14b=zeros(nCellx*nCelly,1);
kk15b=zeros(nCellx*nCelly,1);
DD_N2A=zeros(nCellx*nCelly,1);
DD_NOA=zeros(nCellx*nCelly,1);


% calculate rate and diffusion coefficients at each cell

for i=1:nCellx*nCelly
    
% calculate forward rate coefficients at each cell
kk1f(i,1)=Aaf(1,1)*(T(i,1)^nf(1,1))*exp(-Eaf(1,1)/(k*T(i,1))); 
kk2f(i,1)=Aaf(2,1)*(T(i,1)^nf(2,1))*exp(-Eaf(2,1)/(k*T(i,1))); 
kk3f(i,1)=Aaf(3,1)*(T(i,1)^nf(3,1))*exp(-Eaf(3,1)/(k*T(i,1))); 
kk4f(i,1)=Aaf(4,1)*(T(i,1)^nf(4,1))*exp(-Eaf(4,1)/(k*T(i,1))); 
kk5f(i,1)=Aaf(5,1)*(T(i,1)^nf(5,1))*exp(-Eaf(5,1)/(k*T(i,1))); 
kk6f(i,1)=Aaf(6,1)*(T(i,1)^nf(6,1))*exp(-Eaf(6,1)/(k*T(i,1))); 
kk7f(i,1)=Aaf(7,1)*(T(i,1)^nf(7,1))*exp(-Eaf(7,1)/(k*T(i,1))); 
kk8f(i,1)=Aaf(8,1)*(T(i,1)^nf(8,1))*exp(-Eaf(8,1)/(k*T(i,1))); 
kk9f(i,1)=Aaf(9,1)*(T(i,1)^nf(9,1))*exp(-Eaf(9,1)/(k*T(i,1))); 
kk10f(i,1)=Aaf(10,1)*(T(i,1)^nf(10,1))*exp(-Eaf(10,1)/(k*T(i,1))); 
kk11f(i,1)=Aaf(11,1)*(T(i,1)^nf(11,1))*exp(-Eaf(11,1)/(k*T(i,1)));
kk12f(i,1)=Aaf(12,1)*(T(i,1)^nf(12,1))*exp(-Eaf(12,1)/(k*T(i,1))); 
kk13f(i,1)=Aaf(13,1)*(T(i,1)^nf(13,1))*exp(-Eaf(13,1)/(k*T(i,1))); 
kk14f(i,1)=Aaf(14,1)*(T(i,1)^nf(14,1))*exp(-Eaf(14,1)/(k*T(i,1))); 
kk15f(i,1)=Aaf(15,1)*(T(i,1)^nf(15,1))*exp(-Eaf(15,1)/(k*T(i,1))); 
kk16f(i,1)=Aaf(16,1)*(T(i,1)^nf(16,1))*exp(-Eaf(16,1)/(k*T(i,1)));  
kk17f(i,1)=Aaf(17,1)*(T(i,1)^nf(17,1))*exp(-Eaf(17,1)/(k*T(i,1)));  



% calculate backward rate coefficients at each cell
kk1b(i,1)=Aab(1,1)*(T(i,1)^nb(1,1))*exp(-Eab(1,1)/(k*T(i,1))); 
kk2b(i,1)=Aab(2,1)*(T(i,1)^nb(2,1))*exp(-Eab(2,1)/(k*T(i,1)));
kk3b(i,1)=Aab(3,1)*(T(i,1)^nb(3,1))*exp(-Eab(3,1)/(k*T(i,1))); 
kk4b(i,1)=Aab(4,1)*(T(i,1)^nb(4,1))*exp(-Eab(4,1)/(k*T(i,1))); 
kk5b(i,1)=Aab(5,1)*(T(i,1)^nb(5,1))*exp(-Eab(5,1)/(k*T(i,1)));
kk6b(i,1)=Aab(6,1)*(T(i,1)^nb(6,1))*exp(-Eab(6,1)/(k*T(i,1))); 
kk7b(i,1)=Aab(7,1)*(T(i,1)^nb(7,1))*exp(-Eab(7,1)/(k*T(i,1)));
kk8b(i,1)=Aab(8,1)*(T(i,1)^nb(8,1))*exp(-Eab(8,1)/(k*T(i,1))); 
kk9b(i,1)=Aab(9,1)*(T(i,1)^nb(9,1))*exp(-Eab(9,1)/(k*T(i,1))); 
kk10b(i,1)=Aab(10,1)*(T(i,1)^nb(10,1))*exp(-Eab(10,1)/(k*T(i,1)));  
kk11b(i,1)=Aab(11,1)*(T(i,1)^nb(11,1))*exp(-Eab(11,1)/(k*T(i,1))); 
kk12b(i,1)=Aab(12,1)*(T(i,1)^nb(12,1))*exp(-Eab(12,1)/(k*T(i,1))); 
kk13b(i,1)=Aab(13,1)*(T(i,1)^nb(13,1))*exp(-Eab(13,1)/(k*T(i,1))); 
kk14b(i,1)=Aab(14,1)*(T(i,1)^nb(14,1))*exp(-Eab(14,1)/(k*T(i,1))); 
kk15b(i,1)=Aab(15,1)*(T(i,1)^nb(15,1))*exp(-Eab(15,1)/(k*T(i,1))); 

DD_N2A(i,1)=(1e-4)*(2.63e-3)*sqrt((T(i,1)^3)*((Molar_air+Molar_N2)/(Molar_air*Molar_N2*2)))...
    /(Patm(i,1)*(0.5*(diam_N2+diam_air))^2);

DD_NOA(i,1)=(1e-4)*(2.63e-3)*sqrt((T(i,1)^3)*((Molar_air+Molar_NO)/(Molar_air*Molar_NO*2)))...
    /(Patm(i,1)*(0.5*(diam_NO+diam_air))^2);

    
end

% calculate diffusion coefficients




%% map from DSMC output to Finite Volume Grid
k1f=zeros(nCellx,nCelly);
k2f=zeros(nCellx,nCelly);
k3f=zeros(nCellx,nCelly);
k4f=zeros(nCellx,nCelly);
k5f=zeros(nCellx,nCelly);
k6f=zeros(nCellx,nCelly);
k7f=zeros(nCellx,nCelly);
k8f=zeros(nCellx,nCelly);
k9f=zeros(nCellx,nCelly);
k10f=zeros(nCellx,nCelly);
k11f=zeros(nCellx,nCelly);  
k12f=zeros(nCellx,nCelly);
k13f=zeros(nCellx,nCelly);
k14f=zeros(nCellx,nCelly);
k15f=zeros(nCellx,nCelly);
k16f=zeros(nCellx,nCelly);   
k17f=zeros(nCellx,nCelly); 

k1b=zeros(nCellx,nCelly);
k2b=zeros(nCellx,nCelly);
k3b=zeros(nCellx,nCelly);
k4b=zeros(nCellx,nCelly);
k5b=zeros(nCellx,nCelly);  
k6b=zeros(nCellx,nCelly);
k7b=zeros(nCellx,nCelly);
k8b=zeros(nCellx,nCelly);
k9b=zeros(nCellx,nCelly);
k10b=zeros(nCellx,nCelly);  
k11b=zeros(nCellx,nCelly); 
k12b=zeros(nCellx,nCelly);
k13b=zeros(nCellx,nCelly);
k14b=zeros(nCellx,nCelly);
k15b=zeros(nCellx,nCelly);

for j=2:nCelly+1
    for i=1:nCellx
    coordX(i,j)=data((j-2)*nCellx+i,1); 
    coordY(i,j)=data((j-2)*nCellx+i,2); 
    u(i,j)=data((j-2)*nCellx+i,3); 
    v(i,j)=data((j-2)*nCellx+i,4); 
    rho_N(i,j)=r_N((j-2)*nCellx+i,1);
    rho_O(i,j)=r_O((j-2)*nCellx+i,1);
    rho_N2(i,j)=r_N2((j-2)*nCellx+i,1);
    rho_O2(i,j)=r_O2((j-2)*nCellx+i,1);
    rho_NO(i,j)=r_NO((j-2)*nCellx+i,1);
    
    N_N(i,j)=NN_N((j-2)*nCellx+i,1);
    N_O(i,j)=NN_O((j-2)*nCellx+i,1);
    N_N2(i,j)=NN_N2((j-2)*nCellx+i,1);
    N_O2(i,j)=NN_O2((j-2)*nCellx+i,1);
    N_NO(i,j)=NN_NO((j-2)*nCellx+i,1);
    D_N2A(i,j)=DD_N2A((j-2)*nCellx+i,1);
    D_NOA(i,j)=DD_NOA((j-2)*nCellx+i,1);
    
    
    k1f(i,j)=kk1f((j-2)*nCellx+i,1);
    k2f(i,j)=kk2f((j-2)*nCellx+i,1);
    k3f(i,j)=kk3f((j-2)*nCellx+i,1);
    k4f(i,j)=kk4f((j-2)*nCellx+i,1);
    k5f(i,j)=kk5f((j-2)*nCellx+i,1);
    k6f(i,j)=kk6f((j-2)*nCellx+i,1);
    k7f(i,j)=kk7f((j-2)*nCellx+i,1);
    k8f(i,j)=kk8f((j-2)*nCellx+i,1);
    k9f(i,j)=kk9f((j-2)*nCellx+i,1);
    k10f(i,j)=kk10f((j-2)*nCellx+i,1);
    k11f(i,j)=kk11f((j-2)*nCellx+i,1);  
    k12f(i,j)=kk12f((j-2)*nCellx+i,1);
    k13f(i,j)=kk13f((j-2)*nCellx+i,1);
    k14f(i,j)=kk14f((j-2)*nCellx+i,1);
    k15f(i,j)=kk15f((j-2)*nCellx+i,1);
    k16f(i,j)=kk16f((j-2)*nCellx+i,1);   
    k17f(i,j)=kk17f((j-2)*nCellx+i,1);  
    
    k1b(i,j)=kk1b((j-2)*nCellx+i,1);
    k2b(i,j)=kk2b((j-2)*nCellx+i,1);
    k3b(i,j)=kk3b((j-2)*nCellx+i,1);
    k4b(i,j)=kk4b((j-2)*nCellx+i,1);
    k5b(i,j)=kk5b((j-2)*nCellx+i,1);  
    k6b(i,j)=kk6b((j-2)*nCellx+i,1);
    k7b(i,j)=kk7b((j-2)*nCellx+i,1);
    k8b(i,j)=kk8b((j-2)*nCellx+i,1);
    k9b(i,j)=kk9b((j-2)*nCellx+i,1);
    k10b(i,j)=kk10b((j-2)*nCellx+i,1);  
    k11b(i,j)=kk11b((j-2)*nCellx+i,1); 
    k12b(i,j)=kk12b((j-2)*nCellx+i,1);
    k13b(i,j)=kk13b((j-2)*nCellx+i,1);
    k14b(i,j)=kk14b((j-2)*nCellx+i,1);
    k15b(i,j)=kk15b((j-2)*nCellx+i,1);
    

    end
end

%% Find wall points

% find the wall points for cylinder
wallPoints=zeros(nCellx,nCelly+1);

for i=1:nCellx
    for j=2:nCelly
    
    r=sqrt((-coordX(i,j)+0.015875)*(-coordX(i,j)+0.015875)+coordY(i,j)*coordY(i,j));    
    
    if r<Radius
        wallPoints(i,j)=1;
    else
        wallPoints(i,j)=0;
        
    end      
        
        
    end
end


% Input symmetry cell values
u(:,1)=u(:,2);
v(:,1)=v(:,2);




%% Explicit time loop

for t=1:nTime
    t
    %% inner points
    for i=3:nCellx-1
        for j=3:nCelly
            
           % check for wall points
           
           if wallPoints(i,j)==1
               continue               
           end
           %calculate source terms     
             
           S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms - Second order upwind
                
           conv_N2A(i,j)=-(max(u(i,j),0)/abs(u(i,j)))*(1.5*rho_N2A_old(i,j)*u(i,j)-2*rho_N2A_old(i-1,j)*u(i-1,j)+0.5*rho_N2A_old(i-2,j)*u(i-2,j))/(dx)...
                    -(max(-u(i,j),0)/abs(u(i,j)))*(-0.5*rho_N2A_old(i+2,j)*u(i+2,j)+2*rho_N2A_old(i+1,j)*u(i+1,j)-1.5*rho_N2A_old(i,j)*u(i,j))/(dx)...
                    -(max(v(i,j),0)/abs(v(i,j)))*(1.5*rho_N2A_old(i,j)*v(i,j)-2*rho_N2A_old(i,j-1)*v(i,j-1)+0.5*rho_N2A_old(i,j-2)*v(i,j-2))/(dy)...
                    -(max(-v(i,j),0)/abs(v(i,j)))*(-0.5*rho_N2A_old(i,j+2)*v(i,j+2)+2*rho_N2A_old(i,j+1)*v(i,j+1)-1.5*rho_N2A_old(i,j)*v(i,j))/(dy);     
                
           conv_NOA(i,j)=-(max(u(i,j),0)/abs(u(i,j)))*(1.5*rho_NOA_old(i,j)*u(i,j)-2*rho_NOA_old(i-1,j)*u(i-1,j)+0.5*rho_NOA_old(i-2,j)*u(i-2,j))/(dx)...
                    -(max(-u(i,j),0)/abs(u(i,j)))*(-0.5*rho_NOA_old(i+2,j)*u(i+2,j)+2*rho_NOA_old(i+1,j)*u(i+1,j)-1.5*rho_NOA_old(i,j)*u(i,j))/(dx)...
                    -(max(v(i,j),0)/abs(v(i,j)))*(1.5*rho_NOA_old(i,j)*v(i,j)-2*rho_NOA_old(i,j-1)*v(i,j-1)+0.5*rho_NOA_old(i,j-2)*v(i,j-2))/(dy)...
                    -(max(-v(i,j),0)/abs(v(i,j)))*(-0.5*rho_NOA_old(i,j+2)*v(i,j+2)+2*rho_NOA_old(i,j+1)*v(i,j+1)-1.5*rho_NOA_old(i,j)*v(i,j))/(dy);     
            

                
           % calculate diffusive terms - central differencing
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i-1,j)-2*rho_N2A_old(i,j)+rho_N2A_old(i+1,j))/(dx^2))...
                        +D_N2A(i,j)*((rho_N2A_old(i,j-1)-2*rho_N2A_old(i,j)+rho_N2A_old(i,j+1))/(dy^2)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i-1,j)-2*rho_NOA_old(i,j)+rho_NOA_old(i+1,j))/(dx^2))...
                        +D_NOA(i,j)*((rho_NOA_old(i,j-1)-2*rho_NOA_old(i,j)+rho_NOA_old(i,j+1))/(dy^2)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO;
           
        end
    end
  
    
  %% outlet
  
  for j=3:nCelly
      i=nCellx;
        if wallPoints(i,j)==1
               continue               
        end
           S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms
           
           conv_N2A(i,j)=-(rho_N2A_old(i,j)*u(i,j)-rho_N2A_old(i-1,j)*u(i-1,j))/(dx)...
                    -(rho_N2A_old(i,j)*v(i,j)-rho_N2A_old(i,j-1)*v(i,j-1))/(dy);
                
           conv_NOA(i,j)=-(rho_NOA_old(i,j)*u(i,j)-rho_NOA_old(i,j)*u(i-1,j))/(dx)...
                    -(rho_NOA_old(i,j)*v(i,j)-rho_NOA_old(i,j-1)*v(i,j-1))/(dy);    
                
   
           % calculate diffusive terms
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i,j)-rho_N2A_old(i-1,j))/(dx))...
                    +D_N2A(i,j)*((rho_N2A_old(i,j-1)-2*rho_N2A_old(i,j)+rho_N2A_old(i,j+1))/(dy^2)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i,j)-rho_NOA_old(i-1,j))/(dx))...
                    +D_NOA(i,j)*((rho_NOA_old(i,j-1)-2*rho_NOA_old(i,j)+rho_NOA_old(i,j+1))/(dy^2)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO; 
      
      
  end
%% Free stream
   for i=3:nCellx-1
   j=nCelly+1;            

            %calculate source terms            
             
           S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms
           
           conv_N2A(i,j)=-(rho_N2A_old(i,j)*u(i,j)-rho_N2A_old(i-1,j)*u(i-1,j))/(2*dx)...
                    -(rho_N2A_old(i,j)*v(i,j)-rho_N2A_old(i,j-1)*v(i,j-1))/(dy);
                
           conv_NOA(i,j)=-(rho_NOA_old(i,j)*u(i,j)-rho_NOA_old(i-1,j)*u(i-1,j))/(2*dx)...
                    -(rho_NOA_old(i,j)*v(i,j)-rho_NOA_old(i,j-1)*v(i,j-1))/(dy);    
                
 
           % calculate diffusive terms
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i-1,j)-2*rho_N2A_old(i,j)+rho_N2A_old(i+1,j))/(dx^2))...
                    +D_N2A(i,j)*((rho_N2A_old(i,j)-rho_N2A_old(i,j-1))/(dy)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i-1,j)-2*rho_NOA_old(i,j)+rho_NOA_old(i+1,j))/(dx^2))...
                    +D_NOA(i,j)*((rho_NOA_old(i,j)-rho_NOA_old(i,j-1))/(dy)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO;
   
   
   end
%% Symmetry

for i=3:nCellx-1
    j=2;
       
        if wallPoints(i,j)==1
               continue               
        end
    %calculate source terms

            S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms
           
           conv_N2A(i,j)=-(max(u(i,j),0)/abs(u(i,j)))*(rho_N2A_old(i,j)*u(i,j)-rho_N2A_old(i-1,j)*u(i-1,j))/(dx)...
                    -(max(-u(i,j),0)/abs(u(i,j)))*(rho_N2A_old(i+1,j)*u(i+1,j)-rho_N2A_old(i,j)*u(i,j))/(dx)...
                    -(rho_N2A_old(i,j)*v(i,j)-rho_N2A_old(i,j)*v(i,j))/(dy);
                
           conv_NOA(i,j)=-(max(u(i,j),0)/abs(u(i,j)))*(rho_NOA_old(i,j)*u(i,j)-rho_NOA_old(i-1,j)*u(i-1,j))/(dx)...
                    -(max(-u(i,j),0)/abs(u(i,j)))*(rho_NOA_old(i+1,j)*u(i+1,j)-rho_NOA_old(i,j)*u(i,j))/(dx)...
                    -(rho_NOA_old(i,j)*v(i,j)-rho_NOA_old(i,j)*v(i,j))/(dy);    

           % calculate diffusive terms
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i-1,j)-2*rho_N2A_old(i,j)+rho_N2A_old(i+1,j))/(dx^2))...
                    +D_N2A(i,j)*((rho_N2A_old(i,j)-2*rho_N2A_old(i,j)+rho_N2A_old(i,j+1))/(dy^2)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i-1,j)-2*rho_NOA_old(i,j)+rho_NOA_old(i+1,j))/(dx^2))...
                    +D_NOA(i,j)*((rho_NOA_old(i,j)-2*rho_NOA_old(i,j)+rho_NOA_old(i,j+1))/(dy^2)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO;
    
       
    
end


%% Upper right corner
i=nCellx;
j=nCelly+1;
           %calculate source terms

             
            S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms
           
           conv_N2A(i,j)=-(rho_N2A_old(i,j)*u(i,j)-rho_N2A_old(i-1,j)*u(i-1,j))/(dx)...
                    -(rho_N2A_old(i,j)*v(i,j)-rho_N2A_old(i,j-1)*v(i,j-1))/(dy);
                
           conv_NOA(i,j)=-(rho_NOA_old(i,j)*u(i,j)-rho_NOA_old(i-1,j)*u(i-1,j))/(dx)...
                    -(rho_NOA_old(i,j)*v(i,j)-rho_NOA_old(i,j-1)*v(i,j-1))/(dy);    

%                 
           % calculate diffusive terms
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i,j)-rho_N2A_old(i-1,j))/(dx))...
                    +D_N2A(i,j)*((rho_N2A_old(i,j)-rho_N2A_old(i,j-1))/(dy)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i,j)-rho_NOA_old(i-1,j))/(dx))...
                    +D_NOA(i,j)*((rho_NOA_old(i,j)-rho_NOA_old(i,j-1))/(dy)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO;
%% Lower right corner

i=nCellx;
j=2;


             
           S_N2A(i,j)=k3f(i,j)*N_N(i,j)*N_N(i,j)-k3b(i,j)*N_N2A_old(i,j)-k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)+k5b(i,j)*N_N2(i,j)*N_NOA_old(i,j)...
                      +k4f(i,j)*N_N(i,j)*N_N(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k4b(i,j)*N_N2A_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      +k6f(i,j)*N_N2(i,j)*N_N(i,j)-k6b(i,j)*N_N2A_old(i,j)*N_N(i,j)+k7f(i,j)*N_N2(i,j)*N_N2(i,j)-k7b(i,j)*N_N2(i,j)*N_N2A_old(i,j)...
                      -k9f(i,j)*N_N2A_old(i,j)*N_O2(i,j)+k9b(i,j)*N_N2(i,j)*N_O2(i,j)-k10f(i,j)*N_N2A_old(i,j)*N_O(i,j)+k10b(i,j)*N_N2(i,j)*N_O(i,j)...
                      -k17f(i,j)*N_N2A_old(i,j);
             
           S_NOA(i,j)=k1f(i,j)*N_N(i,j)*N_O(i,j)-k1b(i,j)*N_NOA_old(i,j)+k5f(i,j)*N_N2A_old(i,j)*N_NO(i,j)-k5b(i,j)*N_NOA_old(i,j)*N_N2(i,j)...
                      +k2f(i,j)*N_N(i,j)*N_O(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k2b(i,j)*N_NOA_old(i,j)*(N_N(i,j)+N_O(i,j)+N_N2(i,j)+N_O2(i,j)+N_NO(i,j))...
                      -k11f(i,j)*N_NOA_old(i,j)*N_N2(i,j)+k11b(i,j)*N_N2(i,j)*N_NO(i,j)-k12f(i,j)*N_NOA_old(i,j)*N_O2(i,j)+k12b(i,j)*N_NO(i,j)*N_O2(i,j)...
                      -k13f(i,j)*N_NOA_old(i,j)*N_O(i,j)+k13b(i,j)*N_NO(i,j)*N_O(i,j)-k14f(i,j)*N_NOA_old(i,j)*N_NO(i,j)+k14b(i,j)*N_NO(i,j)*N_NO(i,j)...
                      -k15f(i,j)*N_NOA_old(i,j)*N_N(i,j)+k15b(i,j)*N_NO(i,j)*N_N(i,j)-k16f(i,j)*N_NOA_old(i,j);
             
           % calculate convective terms
           
           conv_N2A(i,j)=-(rho_N2A_old(i,j)*u(i,j)-rho_N2A_old(i-1,j)*u(i-1,j))/(dx)...
                    -(rho_N2A_old(i,j)*v(i,j)-rho_N2A_old(i,j)*v(i,j))/(dy);
                
           conv_NOA(i,j)=-(rho_NOA_old(i,j)*u(i,j)-rho_NOA_old(i-1,j)*u(i-1,j))/(dx)...
                    -(rho_NOA_old(i,j)*v(i,j)-rho_NOA_old(i,j)*v(i,j))/(dy);    
                

           
           % calculate diffusive terms
           diff_N2A(i,j)=D_N2A(i,j)*((rho_N2A_old(i,j)-rho_N2A_old(i-1,j))/(dx))...
                    +D_N2A(i,j)*((rho_N2A_old(i,j)-2*rho_N2A_old(i,j)+rho_N2A_old(i,j+1))/(dy^2)); 
                
           diff_NOA(i,j)=D_NOA(i,j)*((rho_NOA_old(i,j)-rho_NOA_old(i-1,j))/(dx))...
                    +D_NOA(i,j)*((rho_NOA_old(i,j)-2*rho_NOA_old(i,j)+rho_NOA_old(i,j+1))/(dy^2)); 
           
             
           rho_N2A_new(i,j)=rho_N2A_old(i,j)+dt*conv_N2A(i,j)+dt*diff_N2A(i,j)+dt*S_N2A(i,j)*k/R_N2;
           rho_NOA_new(i,j)=rho_NOA_old(i,j)+dt*conv_NOA(i,j)+dt*diff_NOA(i,j)+dt*S_NOA(i,j)*k/R_NO;

           %Fill in the symmetry cells
           rho_N2A_new(:,1)=rho_N2A_new(:,2);
           rho_NOA_new(:,1)=rho_NOA_new(:,2);  
    

      
  %% calculate errors and update the solution for densities of N2A and NOA 
   
   rho_N2A_new(isnan(rho_N2A_new))=0;
   rho_NOA_new(isnan(rho_NOA_new))=0;
  err_N2A=abs(rho_N2A_old-rho_N2A_new)./rho_N2A_old;
   err_NOA=abs(rho_NOA_old-rho_NOA_new)./rho_NOA_old;
    
   
   
   rho_N2A_old=rho_N2A_new;
   rho_NOA_old=rho_NOA_new;
   N_N2A_old=rho_N2A_old*R_N2/k;
   N_NOA_old=rho_NOA_old*R_NO/k;
   
   
   NOA70(t,1)=N_NOA_old(140,2);
   N2A70(t,1)=N_N2A_old(140,2);
%    SNOA70(t,1)=S_NOA(70,2)*k/R_NO;
%    SN2A70(t,1)=S_N2A(70,2)*k/R_N2;
%    convNOA70(t,1)=conv_NOA(70,2);
%    convN2A70(t,1)=conv_N2A(70,2);
%    diffNOA70(t,1)=diff_NOA(70,2);
%    diffN2A70(t,1)=diff_N2A(70,2);
    
end
N2A=zeros(nCellx*nCelly,1);
NOA=zeros(nCellx*nCelly,1);
%% Convert to vector forms and write output for tecplot            
for j=2:nCelly+1
    for i=1:nCellx
    
 N2A((j-2)*nCellx+i,1)=N_N2A_old(i,j); 
 NOA((j-2)*nCellx+i,1)=N_NOA_old(i,j);
 XX((j-2)*nCellx+i,1)=coordX(i,j);
 YY((j-2)*nCellx+i,1)=coordY(i,j);
 uu((j-2)*nCellx+i,1)=u(i,j);
 VV((j-2)*nCellx+i,1)=v(i,j);  
    end
end            
time=dt:dt:dt*nTime;
time=time';
N2Acm3=N2A*1e-6;
NOAcm3=NOA*1e-6;

timeOutput=[time,NOA70,N2A70];
% timeOutputTerms=[time,SNOA70,SN2A70,convNOA70,convN2A70,diffNOA70,diffN2A70];
Outputm3=[XX,YY,uu,VV,NOA,N2A];
% Outputcm3=[XX,YY,NOAcm3,N2Acm3];
save MassTRFullSetm3.dat Outputm3 -ascii   
% save MassTRNOAN2Adensitycm3.dat Outputcm3 -ascii   
save MassTRNOAN2Adensitytime.dat timeOutput -ascii 
% save Termstime.dat timeOutputTerms -ascii 
            
            
            
            
            
            
            
            
            




