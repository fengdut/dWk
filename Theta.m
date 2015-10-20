function [Theta_3D] = Theta(b_lambda_3D,ntheta,dtheta,nx,xarray,nL,Larray)
global eps
Theta_3D(1:ntheta,1:nx,1:nL) =0;

for iL=1:nL
    for ix=1:nx
        for itheta=2:ntheta
           Theta_3D(itheta,ix,iL) =  Theta_3D(itheta-1,ix,iL) +b_lambda_3D(itheta,ix,iL)*dtheta*0.5 +b_lambda_3D(itheta-1,ix,iL)*dtheta*0.5 ;
        end
    end
end


% for iL=1:nL
%     parfor ix=1:nx
%         for itheta=2:ntheta
%             tc=eps*xarray(ix);
%            fun=@(t)1./((1-tc*cos(t)).*sqrt(1-Larray(iL)./(1-tc*cos(t))));
%            Theta_3D(itheta,ix,iL) = quad(fun,(itheta-2)*dtheta,(itheta-1)*dtheta);
%         end
%     end
% end

 for iL=1:nL
     for ix=1:nx
        for itheta=3:ntheta
           Theta_3D(itheta,ix,iL) =  Theta_3D(itheta-1,ix,iL)+ Theta_3D(itheta,ix,iL);
        end
    end
 end

end