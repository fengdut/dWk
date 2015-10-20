function [F_3D,F_x,F_L,F_E]=F0_3D(nx,xarray,nL,Larray,nE,Earray,x0,xd,L0,Ld,E0,Ed,Ec)
%slowing down distribution function
F_3D(1:nx,1:nL,1:nE)=0;

for iE=1:nE
    for iL=1:nL
        for ix=1:nx
            tE=Earray(iE);
            tL=Larray(iL);
            tx=xarray(ix);
            ef=1/(tE^1.5 +Ec^1.5) *erfc((tE-E0)/Ed);
            Lf=exp(-((tL-L0)/Ld)^2);
            xf=exp(-((tx-x0)/xd)^2);
            F_3D(ix,iL,iE) = ef*Lf*xf;
        end
    end
end



end