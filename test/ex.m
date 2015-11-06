clear
close all
x=1e-06;
Lambda=1e-06;
t=1e-6;
%b(t,x)*sqrt(1-Lambda./b(t,x))
kappa=(1-Lambda.*(1-eps1(x)))./(2.0.*eps1(x).*Lambda);
chi=pi*sqrt(kappa.*eps1(x).*Lambda/2)./...
    ellipke(1./kappa)./...
    (b(t,x)*sqrt(1-Lambda./b(t,x)));
Gfun(Lambda,t,x)
nt=10;
dt=2*pi/nt;
t=0:dt:2*pi;
p=1;
fi=Y1(Lambda,t,t,x,p);
simp(length(t),dt,fi)
