%this is subroutine to calculate lambda3d(1:ntheta,1:nx,1:nL)

function lambda3d=lambda(ntheta,thetaarray,nx,xarray,nL,Larray)
lambda3d(1:ntheta,1:nx,1:nL)=0;

for ix=1:nx
    for iL=1:nL
        for itheta=1:ntheta
            bt=b(thetaarray(itheta),xarray(ix));
            lambda3d(itheta,ix,iL) = 1/(tb*sqrt(1-))
        end
    end
end
        

end