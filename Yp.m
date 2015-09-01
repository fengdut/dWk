function f=Yp(Lambda,epsilon,x,tau,p) 

% 
% f=quad(@(t)Y1(Lambda,epsilon,t,x,tau,p),0,2*pi-0.1)*...
%    omega_b(Lambda,epsilon,x)/2/pi;

f=quad(@(t)Y1(Lambda,epsilon,t,x,tau,p),0,2*pi-0.1);
 
f=f*    omega_b(Lambda,epsilon,x)/2/pi;


end
