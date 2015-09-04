function f=tfun(Lambda,t,x,tau) % omega_b*t 
% only one of Lambda,t,x allowed to be vector, two vectors present same time not permitted for muti-variables function
kappa=(1-Lambda.*(1-eps1(x)))./(2.0.*eps1(x).*Lambda);
ta=0.0;
nt=10;% even number gives better accuracy for Simpson method
dt=(t-ta)/nt;
tt=ta:dt:t;
fi=dtfun(Lambda,tt,x);
f=tau*pi*sqrt(kappa.*eps1(x).*Lambda/2).*...
    simp(length(tt),dt,fi)./...
    ellipke(1./kappa);
end
