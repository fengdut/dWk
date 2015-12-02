clear;
%%%%%%%%%%%%%%%% tokamak parameters %%%%%%%%%%%%%
global R0 a eps
R0=1.65;
a=0.4;
eps=a/R0;

%%%%%%%%%%%%%%%% equilibrum parameters %%%%%%%%%%
global rs
rs=0.5;    %%should be consistent with q profile

%%%%%%%%%%%%%% mode parameters %%%%%%%%%%%%%%%%% 
global n_mode m_mode m_p delta_r
n_mode=1;  %toroidal mode number
m_mode=1;  %poloidal mode number
m_p=0;
delta_r=0.06;

%%%%%%%%%%%%% calculation mesh and boundary %%%%%%
global nx nL nE ntheta xa xb La Lb Ea Eb dx dL dE dtheta 
global xarray Larray Earray thetaarray
nx=100;        %grid number of radial
nL=100;        %grid number of lambda 
nE=100;        %grid number of energy
ntheta=61;     %grid number of theta

xa=1e-6;     %r_min, avoid r==0
xb=1.0;      %r_max
Ea=0.1;      %E_min
Eb=5.0;      %E_max
La=1e-6;      %Lambda_min
Lb=0.1;       %lambda_max

dx=(xb-xa)/(nx-1);
dL=(Lb-La)/(nL-1);
dE=(Eb-Ea)/(nE-1);
dtheta=2*pi/(ntheta-1);
xarray=xa:dx:xb;
Larray=La:dL:Lb;
Earray=Ea:dE:Eb;
thetaarray=0:dtheta:2*pi;

if(mod(nx,3)~=1)
    display('grid number must be 3*n+1');
    display('exit the code');
    return;
end

%%%%%%%energetic particles parameters%%%%%%
global rd r0 Ld L0 Ed E0 Ec rhoh sigma 
rd=0.3;         %Delta r
Ld=0.02;        %Delta Lambda
r0=0;           %r_0
L0=0.05;        %Lambda_0
E0=Eb;          %E_0
Ed=1;           %Delta E
Ec=0.8*Eb;      %E_c
rhoh=0.08;
sigma=1;

%%%%%%%%% calcuate non-omega parts %%%%%
Non_omega_parts;             


%%%%%%%%% 
tic
display('begin to calculate delta wk');
delta_wk=dwk(1.5+0.001i);
toc
delta_wk



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

