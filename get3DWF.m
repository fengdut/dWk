
tic;

small_imag=0.1e-2i;


n_omega=1;
nn_omega=100;

domega=0.02;
omega_0=1.5;

omega_a=0;
sum_dwk=0;
omega_a(1:n_omega)=0+0i;
sum_dwk(1:n_omega)=0;

dsmall=0.001i;
for iiomega=1:nn_omega
%small_imag=small_imag +dsmall;
for iomega=1:n_omega

    omega=omega_0+domega*(iomega-1)+small_imag*iiomega;
omega_a(iiomega)=omega;
m_p=0;

R_3D=n_mode *omega_phi_3D +m_p*omega_b_3D -omega; %R_3D(r,Lambda,E)             %R_3D(x,L,E)

Yp_R_3D(1:nx,1:nL,1:nE) = 0;                                                %Yp_R_3D(x,L,E)
for ix=1:nx
    for iL=1:nL
        for iE=1:nE
            Yp_R_3D(ix,iL,iE)=Yp2(ix,iL)/R_3D(ix,iL,iE);
        end
    end
end

WF3D(1:nx,1:nL,1:nE)=0+0i;                                                     %WF3D(x,L,E)
for ix=1:nx
    for iL=1:nL
        for iE=1:nE
           WF3D(ix,iL,iE) =J_q(ix)*Earray(iE)*tau_b_3D(ix,iL,iE)*(F_E_3D(ix,iL,iE)*real(omega)-dFdE_omega_star_3D(ix,iL,iE))*Yp_R_3D(ix,iL,iE);
% WF3D(ix,iL,iE) =Yp_R_3D(ix,iL,iE);
        end
    end
end

sum_dwk(iiomega)=simpintegral_3D(WF3D,nx,dx,nL,dL,nE,dE);

end



end
xp=imag(omega_a);
yp=imag(sum_dwk);
plot(imag(omega_a),imag(sum_dwk),'o--');
hold all;
toc








