function F_E_3D=dFdE(nx,xarray,nL,Larray,nE,Earray) % the partial derivative of the distribution 
global rd r0 Ld L0 Ed E0 Ec;
F_E_3D(1:nx,1:nL,1:nE) = 0;


expr=exp(-((xarray-r0)/rd).^2);
expL=exp(-((Larray-L0)/Ld).^2);
expE=exp(-((Earray-E0)/Ed).^2);

E32=Earray.^1.5 +Ec^3/2;
expE= 2*expE./(pi.^0.5*Ed*E32);
erfcE=3*Earray.^0.5 .*erfc((Earray-E0)/Ed)/(2*E32.^2);
exp_erfc_E=expE+erfcE;

for ix=1:nx
    for iL=1:nL
        for iE=1:nE
            F_E_3D(ix,iL,iE) = -exp_erfc_E(iE) *expr(ix) *expL(iL); 
        end
    end
end


end