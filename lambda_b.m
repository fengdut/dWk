function [lambda_b_3D] = lambda_b(ntheta,thetaarray,nx,xarray,nL,Larray)

lambda_b_3D(1:ntheta,1:nx,1:nL) =0;

for iL=1:nL
    for ix=1:nx
        for itheta=1:ntheta
            lambda_b_3D(itheta,ix,iL) = Larray(iL) / (b(thetaarray(itheta),xarray(ix)));
        end
    end
end

end
    