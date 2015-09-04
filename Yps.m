function Yps_2D =Yps(G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,nx,nL,ntheta,dtheta,sigma,p)
Yps_2D(1:nx,1:nL) = 0;

for iL=1:nL
    for ix=1:nx
        sum =0;
        tY(1:ntheta)=0;
        
        for it=1:ntheta
            texp=exp(-1i*Chi_2D(ix,iL)*p*Theta_3D(it,ix,iL));
     %  texp=1;
            tY(it)=G_2D(it,ix)*b_lambda_3D(it,ix,iL)*(lambda_b_3D(it,ix,iL)+2*(1-lambda_b_3D(it,ix,iL)))*texp;
        end
        
        Yps_2D(ix,iL) =simpintegral(tY,ntheta,dtheta)*Chi_2D(ix,iL)/(2*pi);
    end
end

end