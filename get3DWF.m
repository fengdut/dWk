
clear;
format long;
global R0 a rhoh xa xb La Lb Ea Eb dx dL dE
global nx nL nE eps rs
global tfun3d

tfun3d=0;
R0=1.65;
a=0.4;
rhoh=0.05;
rs=0.5;
eps=a/1.65;
delta_r=0.1;


xa=1e-6;
xb=1.0;

Ea=2.0;
Eb=4.0;

La=1e-6;
Lb=0.7;

nx=202;
nL=10;
nE=31; 

ntheta=61;


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
p=-1;

tic

[lambda_b_3D]=lambda_b(ntheta,thetaarray,nx,xarray,nL,Larray);
[b_lambda_3D]=b_lambda(lambda_b_3D,ntheta,thetaarray,nx,xarray,nL,Larray);
[Theta_3D]   =Theta(b_lambda_3D,ntheta,dtheta,nx,xarray,nL,Larray);
[G_2D]       =G(ntheta,thetaarray,nx,xarray,rs,delta_r);

[Chi_2D]     = Chi(nx,xarray,nL,Larray,sigma);
[Yps_2D]     = Yps (G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,nx,nL,ntheta,dtheta,sigma,p);
toc


Yps_2D(1,2)


%Theta_3D(ntheta,1,1)




omega=0.1+0.01i;

WF3D(1:nx,1:nL,1:nE)=0;









