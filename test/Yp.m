function f=Yp(Lambda,tf,x,p) 
global nt
kappa=(1-Lambda.*(1-eps1(x)))./(2.0.*eps1(x).*Lambda);
nt=500; % same as main/test subroutine
T=2*pi;
dt=(T-1e-6)/(nt);
t=1e-6:dt:T;
tf1=zeros(1,length(t));
tf1(:)=tf;
fi=Y1(Lambda,t,tf1,x,p); % tf is array which is same as t array
%f=pi*sqrt(kappa.*eps1(x).*Lambda/2).*...
%    simp(length(t),dt,fi)./...
%    ellipke(1./kappa)./(2*pi);
f=simp(length(t),dt,fi);
end
