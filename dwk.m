function delta_wk=dwk(omega)
global nx nL nE Earray dx dL dE
global n_mode m_p omega_phi_3D omega_b_3D Yp2 J_q
global tau_b_3D F_E_3D dFdE_omega_star_3D

iomega=imag(omega);
%romega=real(omega);
if(abs(iomega)<0.001)
    disploy('small image part of omega');
end
R_3D =n_mode *omega_phi_3D + m_p*omega_b_3D -omega;
Yp_R_3D(1:nx,1:nL,1:nE)=0;

for ix=1:nx
    for iL=1:nL
        for iE=1:nE
            Yp_R_3D(ix,iL,iE)=Yp2(ix,iL)/R_3D(ix,iL,iE);
        end
    end
end

W3D(1:nx,1:nL,1:nE)=0+0i;

for ix=1:nx
    for iL=1:nL
        for iE=1:nE
           W3D(ix,iL,iE) =J_q(ix)*Earray(iE)*tau_b_3D(ix,iL,iE)*(F_E_3D(ix,iL,iE)*real(omega)-dFdE_omega_star_3D(ix,iL,iE))*Yp_R_3D(ix,iL,iE);
            % W3D(ix,iL,iE) =Yp_R_3D(ix,iL,iE);
        end
    end
end

delta_wk = simpintegral_3D(W3D,nx,dx,nL,dL,nE,dE);
end