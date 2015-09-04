function f=Ypi(Lambda,epsilon,x,tau,p) 

% 
% f=quad(@(t)Y1(Lambda,epsilon,t,x,tau,p),0,2*pi-0.1)*...
%    omega_b(Lambda,epsilon,x)/2/pi;
global tomega_b

ntheta=100;
dtheta=2*pi/(ntheta-1);

sum=0;
 for ti=1:ntheta
     ttheta=(ti-1)*dtheta;
     tYp(ti) = Y1(Lambda,epsilon,ttheta,x,tau,p);
     sum=sum+tYp(ti);
 end
 
 for ti=2:ntheta-1
     sum=sum+2*tYp(ti);
 end
 
 for ti=4:3:ntheta-1
     sum = sum-tYp(ti);
 end
 
 f=sum*dtheta*3/8;
 
 
 
%f=quad(@(t)Y1(Lambda,epsilon,t,x,tau,p),0,2*pi,0.01);
 
f=f* tomega_b/(2*pi);


end
