

function F_r_3D=dFdr(nx,xarray,nL,Larray,nE,Earray) % the partial derivative of the distribution 
global rd r0 Ld L0 Ed E0 Ec;
F_r_3D(1:nx,1:nL,1:nE) = 0;


expr=exp(-((xarray-r0)/rd).^2);
expL=exp(-((Larray-L0)/Ld).^2);


E32=Earray.^1.5 +Ec^1.5;
erfcE=erfc((Earray-E0)/Ed)./E32;


tr=(2.*(r0-xarray)/rd^2);

expr=tr.*expr;
for ix=1:nx
    for iL=1:nL
        for iE=1:nE
            F_r_3D(ix,iL,iE) = erfcE(iE) *expr(ix) *expL(iL); 
        end
    end
end


end