function f=WF(Yp,R_term,Jac,qprofile,epsilon,density,dFdE,FL,omega_b,omega_minus_omega_star) % the integrand of delta_W, omega is imaginary number normalized to vth/R0. vth=sqrt(2Th/M)

g=sum(abs(Yp).^2./R_term);
f=Jac./qprofile.*epsilon.^3.*density*...
    dFdE.*FL.*....
  2*pi./omega_b.*...
  omega_minus_omega_star.*g;    
end
