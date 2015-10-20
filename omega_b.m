function omega_b_3D=omega_b(nx,xarray,nE,Earray,nL,Larray,kappa,K,q_1D)

global  eps



omega_b_3D(1:nx,1:nL,1:nE)=0;
for ix=1:nx
    for iL=1:nL
        for iE=1:nE
            q=q_1D(ix);
            omega_b_3D(ix,iL,iE) = pi *sqrt(kappa(ix,iL))/K(ix,iL) * sqrt(xarray(ix)*eps*Larray(iL)/2) /q*sqrt(Earray(iE));
        end
    end
end

end
