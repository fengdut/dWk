function f=omega_b(kappa,qprofile,Lambda,epsilon,eps1,K)
f=pi*sqrt(kappa.*eps1.*epsilon.*Lambda/2.0)./...
    (qprofile.*K); % transit frequency
end
