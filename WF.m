function f=WF(omega,Lambda,epsilon,x,tau) % the integrand of delta_W, omega is imaginary number normalized to vth/R0. vth=sqrt(2Th/M)
g=0.0;
n=1;

if(Lambda>=1-eps1(x)-0.01)
    f =0;
else
for p=-1:1
    g=g+abs(Yp(Lambda,epsilon,x,tau,p))^2/...
    (n*omega_phi(Lambda,epsilon,x)+...
     p*omega_b(Lambda,epsilon,x)-omega);
end
f=Jac(x)/qprofile(x)*epsilon^3*density(x)*...
    pFpE(Lambda,epsilon)*...
  2*pi/omega_b(Lambda,epsilon,x)*...
  (omega-omega_star(epsilon,x))*g;    
end
end
