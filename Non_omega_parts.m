
format long;
global R0 a rhoh xa xb La Lb Ea Eb dx dL dE
global nx nL nE eps rs
global n_mode m_mode m_p

n_mode=1;  %toroidal mode number
m_mode=1;  %poloidal mode number
m_p=0;

R0=1.65;
a=0.4;
rhoh=0.08;
rs=0.5;
eps=a/1.65;
delta_r=0.06;


xa=1e-6;     %r_min, avoid r==0
xb=1.0;      %r_max

Ea=0.1;      %E_min
Eb=5.0;      %E_max

La=1e-6;      %Lambda_min
Lb=0.1;       %lambda_max

nx=100;        %grid number of radial
nL=100;        %grid number of lambda 
nE=100;        %grid number of energy

ntheta=22;     %grid number of theta

modn=mod(nx,3);
if(modn~=1)
    display('grid number must be 3*n+1');
    display('exit the code');
    return;
end

%distribution function
global rd r0 Ld L0 Ed E0 Ec 
global xarray Larray Earray
rd=0.3;         %Delta r
Ld=0.02;        %Delta Lambda
r0=0;           %r_0
L0=0.05;        %Lambda_0
E0=Eb;          %E_0
Ed=1;           %Delta E
Ec=0.8*Eb;      %E_c

dx=(xb-xa)/(nx-1);
dL=(Lb-La)/(nL-1);
dE=(Eb-Ea)/(nE-1);
dtheta=2*pi/(ntheta-1);

thetaarray=0:dtheta:2*pi;
xarray=xa:dx:xb;
Larray=La:dL:Lb;
Earray=Ea:dE:Eb;


tau=1.0;
sigma=1;

tic
global omega_phi_3D  omega_b_3D Yp2 J_q tau_b_3D F_E_3D dFdE_omega_star_3D
q_1D=qprofile(nx,xarray); 
F_3D=F0_3D(nx,xarray,nL,Larray,nE,Earray,r0,rd,L0,Ld,E0,Ed,Ec);
F_E_3D=dFdE(nx,xarray,nL,Larray,nE,Earray);                                 %df_de (x,L,E)
F_r_3D=dFdr(nx,xarray,nL,Larray,nE,Earray);                                 %df_dr(x,L,E)

[lambda_b_3D]=lambda_b(ntheta,thetaarray,nx,xarray,nL,Larray);              %lambda_b_3D(theta,x,L)
toc
tic
[b_lambda_3D]=b_lambda(lambda_b_3D,ntheta,thetaarray,nx,xarray,nL,Larray);  %b_lambda_3D(theta,x,L)
toc
tic
[Theta_3D]   =Theta(b_lambda_3D,ntheta,dtheta,nx,xarray,nL,Larray);         %Theta_3D(theta,x,L)

[G_2D]       =G(ntheta,thetaarray,nx,xarray,rs,delta_r);                    %G_2D(theta,x)
 
[Chi_2D,kappa_2D,K_2D] = Chi(nx,xarray,nL,Larray,sigma);                    %Chi_2D(x,L),kappa_2D(x,L),K_2D(x,L)

[Yps_2D]     = Yps (G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,nx,nL,ntheta,dtheta,sigma,m_p);  %Yps_2D(x,L)


omega_b_3D=omega_b(nx,xarray,nE,Earray,nL,Larray,kappa_2D,K_2D,q_1D);            %omega_b_3D(x,L,E)

omega_phi_3D=omega_phi(omega_b_3D,nx,xarray,q_1D);                               %omega_phi_3D(x,L,E)

Yp2 = abs(Yps_2D).^2; %Yp2(r,lambda)                                             %Yp2(x,L)


tau_b_3D=2*pi./omega_b_3D;                                                  
J_1D=Jac(nx,xarray);

J_q=J_1D./q_1D;
dFdE_omega_star_3D= dFdE_omega_star(nx,xarray,nL,Larray,nE,Earray,F_r_3D);
toc

display('non omega parts finish');




