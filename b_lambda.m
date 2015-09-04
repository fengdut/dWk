function [b_lambda_3D] = b_lambda(lambda_b_3D,ntheta,thetaarray,nx,xarray,nL,Larray)

b_lambda_3D(1:ntheta,1:nx,1:nL) =0;

for iL=1:nL
    for ix=1:nx
        for itheta=1:ntheta
            tb=b(thetaarray(itheta),xarray(ix));
            b_lambda_3D(itheta,ix,iL) = 1 / (tb*sqrt(1-lambda_b_3D(itheta,ix,iL)));
        end
    end
end

end
    