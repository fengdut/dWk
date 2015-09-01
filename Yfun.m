function f=Yfun(Lambda,epsilon,x,tau,p)
rho=0.1;
f=omega_b(Lambda,epsilon,x)/2/pi.*...
    quad(@(t)Y1(Lambda,epsilon,t,x+rho*cos(x),tau,p),0,2*pi);
end
