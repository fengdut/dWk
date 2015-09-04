function [Chi_2D] = Chi(nx,xarray,nL,Larray,sigma)
Chi_2D(1:nx,1:nL)=0;


for iL=1:nL
    for ix=1:nx
        teps=eps1(xarray(ix));
        k=(1-Larray(iL)*(1-teps));
        k=k/(2*Larray(iL)*teps);
        if(k>1)  % only for passing particles
        K=ellipke(1/k);
        Chi_2D(ix,iL) = sigma*pi *sqrt(k *teps *Larray(iL)/2)/K;
        end
    end
end
end