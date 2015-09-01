function f=omega_satr(epsilon,x)
global R0 a rhoh
m=1.0; % poloidal mode number
f=m*dFdE(epsilon).*density(x)./denp(x)/2.0.*R0/a*rhoh/a;
end
