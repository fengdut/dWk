function f=tfun(kappa,Lambda,t,x,eps1,K) % omega_b*t 
% only one of Lambda,t,x allowed to be vector, 
%two vectors present same time not permitted 
%for muti-variables function in general
%kappa=(1-Lambda.*(1-eps1(x)))./(2.0.*eps1(x).*Lambda);
global tau
ta=0.0;
nt1=10;% even number gives better accuracy for Simpson method
dt=(t-ta)/nt1;
tt=ta:dt:t;
%fi=1./(b(t,x).*sqrt(1-Lambda./b(t,x)));
fi=dtfun(Lambda,tt,x);
f=tau*pi*sqrt(kappa.*eps1.*Lambda/2).*...
    simp(length(tt),dt,fi)./K;
end
