global nx nL nE rs
global m_p
global xarray Larray Earray thetaarray
%distribution function
global rd r0 Ld L0 Ed E0 Ec sigma 


tic
global omega_phi_3D  omega_b_3D Yp2 J_q tau_b_3D F_E_3D dFdE_omega_star_3D
q_1D=qprofile(nx,xarray); 
F_3D=F0_3D(nx,xarray,nL,Larray,nE,Earray,r0,rd,L0,Ld,E0,Ed,Ec);
F_E_3D=dFdE(nx,xarray,nL,Larray,nE,Earray);                                 %df_de (x,L,E)
F_r_3D=dFdr(nx,xarray,nL,Larray,nE,Earray);                                 %df_dr(x,L,E)

lambda_b_3D=lambda_b(ntheta,thetaarray,nx,xarray,nL,Larray);              %lambda_b_3D(theta,x,L)
toc
tic
b_lambda_3D=b_lambda(lambda_b_3D,ntheta,thetaarray,nx,xarray,nL,Larray);  %b_lambda_3D(theta,x,L)
toc
tic
Theta_3D   =Theta(b_lambda_3D,ntheta,dtheta,nx,xarray,nL,Larray);         %Theta_3D(theta,x,L)

G_2D       =G(ntheta,thetaarray,nx,xarray,rs,delta_r);                    %G_2D(theta,x)
 
[Chi_2D,kappa_2D,K_2D] = Chi(nx,xarray,nL,Larray,sigma);                    %Chi_2D(x,L),kappa_2D(x,L),K_2D(x,L)

Yps_2D     = Yps (G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,nx,nL,ntheta,dtheta,sigma,m_p);  %Yps_2D(x,L)

omega_b_3D=omega_b(nx,xarray,nE,Earray,nL,Larray,kappa_2D,K_2D,q_1D);            %omega_b_3D(x,L,E)

omega_phi_3D=omega_phi(omega_b_3D,nx,xarray,q_1D);                               %omega_phi_3D(x,L,E)

Yp2 = abs(Yps_2D).^2; %Yp2(r,lambda)                                             %Yp2(x,L)

tau_b_3D=2*pi./omega_b_3D;                                                  
J_1D=Jac(nx,xarray);

J_q=J_1D./q_1D;
dFdE_omega_star_3D= dFdE_omega_star(nx,xarray,nL,Larray,nE,Earray,F_r_3D);
toc

display('non-omega parts finish');




