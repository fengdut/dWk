function omega_phi_3D=omega_phi(omega_b,nx,xarray,q_1D)

%omega_phi_3D(1:nx,1:nL,1:nE)=0;
for ix =1:nx
    q=q_1D(ix);
    omega_phi_3D(ix,:,:)=q*omega_b(ix,:,:);
end
end
