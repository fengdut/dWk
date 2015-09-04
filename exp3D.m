%this is subroutine to calculate exptheta(1:ntheta,1:nx,1:nL)

function exp3Darrry=exp3D(ntheta,thetaarray,nx,xarray,nL,Larray,sigma,p,epsarray)

    exp3Darray(1:ntheta,1:nx,1:nL)=0;
    dtheta=(thetaarray(ntheta)-thetaarray(1))/size(thetaarray,1);
    
    for ix=1:nx
        for iL=1:nL
            for it=1:ntheta
                tb=b(thetaarray(it)),xarray(ix);
                tsum = 1/(tb*(1-Larray(iL)/tb)^0.5);
                exp3Darray(ix,iL,it) =exp3Darray(ix,iL,it)+ tsum;
            end
            end
    end
    
    
    for ix=1:nx
        for iL=1:nL
            for it=1:ntheta
                tb=b(thetaarray(it)),xarray(ix);
                tsum = 1/(tb*(1-Larray(iL)/tb)^0.5);
                exp3Darray(ix,iL,it) =exp3Darray(ix,iL,it)+ tsum;
            end
            end
    end
    
    
    
end


