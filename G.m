function [G_2D] =G(ntheta,tarray,nx,xarray,r_s,delta_r)
global eps;
G_2D(1:ntheta,1:nx)=0;
cost = cos(tarray);
sint = sin(tarray);

xi_r(1:nx) =1;
xi_t(1:nx) =-1;

nx1=1;
nx2=1;
for ix=1:nx
    if(xarray(ix)>=r_s-delta_r/2)
        nx1=ix;
            break;
    end

end

for ix=1:nx
    if(xarray(ix)>=r_s+delta_r/2)
        nx2=ix;
            break;
    end

end


xi_r(nx2:nx)=0;
xi_t(nx2:nx)=0;
for ix=nx1:nx2-1
    xi_r(ix) = (delta_r - xarray(ix)+r_s-delta_r/2)/delta_r;
    xi_t(ix) =-(delta_r - 2*xarray(ix)+r_s-delta_r/2)/delta_r;
end
xi_t =xi_t.*xarray;
xi_r =cos(tarray)'*xi_r;
xi_t =sin(tarray)'*xi_t;

grr(1:ntheta,1:nx)=0;
gtt(1:ntheta,1:nx)=0;
grt(1:ntheta,1:nx)=0;
kappa_r(1:ntheta,1:nx)=0;
kappa_t(1:ntheta,1:nx)=0;
for ix=1:nx
    for it=1:ntheta
        teps=eps1(xarray(ix));
        grr(it,ix) = 1 + teps*cost(it)/2;
        gtt(it,ix) = 1/xarray(ix)^2 *(1 - 5/2 *teps *cost(it));
        grt(it,ix) = -3/(2*xarray(ix))*teps*sint(it);
        kappa_r(it,ix) = -1/(1/eps +xarray(ix)*cost(it))*cost(it);
        kappa_t(it,ix) = teps* sint(it);
    end
end
G_2D = (gtt.*kappa_t +grt.*kappa_r).*xi_t +(grr.*kappa_r +grt.*kappa_t).*xi_r;

end