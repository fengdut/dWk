function f=WF(omega,Lambda,epsilon,x,Yp) % the integrand of delta_W, omega is imaginary number normalized to vth/R0. vth=sqrt(2Th/M)
global nn
g=0.0;
p=-1:1;
np=length(p);
Yp1=zeros(1,np);
Yp1(:)=Yp;
    for j=1:np
       g=g+abs(Yp1(j))^2/...
       (nn*omega_phi(Lambda,epsilon,x)+...
        p(j)*omega_b(Lambda,epsilon,x)-omega);
    end
f=Jac(x)./qprofile(x).*epsilon.^3.*density(x)*...
    dFdE(epsilon).*FL(Lambda).*....
  2*pi./omega_b(Lambda,epsilon,x).*...
  (omega-omega_star(epsilon,x)).*g;    
end
