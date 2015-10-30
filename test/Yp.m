function f=Yp(Lambda,epsilon,tf,x,p,tau) 
global  dt nt ta tb
kappa=(1-Lambda.*(1-eps1(x)))./(2.0.*eps1(x).*Lambda);
tf1=zeros(1,nt+1);
tf1(:)=tf;               % tf is time array as function of theta
t=ta:dt:tb;
fi=Y1(Lambda,epsilon,t,tf1,x,p); 
f=tau*pi*sqrt(kappa.*eps1(x).*Lambda/2).*...
    simp(nt+1,dt,fi)./...
    ellipke(1./kappa)./(2*pi);
%f=simp(length(t),dt,fi);
end
