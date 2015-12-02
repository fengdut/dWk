function [Chi_2D,kappa_2D,K_2D] = Chi(nx,xarray,nL,Larray,sigma)
global eps
Chi_2D(1:nx,1:nL)=0;

kappa_2D(1:nx,1:nL)=0;
K_2D(1:nx,1:nL)=0;
for iL=1:nL
    for ix=1:nx
        teps=eps*(xarray(ix));
        k=(1-Larray(iL)*(1-teps));
        k=k/(2*Larray(iL)*teps);
        kappa_2D(ix,iL)=k;
        if(k>1)  % only for passing particles
        K=ellipke(1/k);
        K_2D(ix,iL)=K;
        Chi_2D(ix,iL) = sigma*pi *sqrt(k *teps *Larray(iL)/2)/K;
        end
    end
end
end