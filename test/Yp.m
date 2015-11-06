function f=Yp(Y1i,kappa,Lambda,eps1,K) 
global  dt nt  tau
f=tau*pi*sqrt(kappa.*eps1.*Lambda/2).*...
    simp(nt+1,dt,Y1i)./...
    K./(2*pi);
%f=simp(length(t),dt,fi);
end
