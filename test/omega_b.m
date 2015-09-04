function f=omega_b(Lambda,epsilon,x)
kappa=(1.0-Lambda.*(1-eps1(x)))./(2.0*eps1(x).*Lambda);
f=pi*sqrt(kappa).*sqrt(eps1(x).*epsilon.*Lambda/2.0)./...
    (qprofile(x).*ellipke(1./kappa)); % transit frequency
end
