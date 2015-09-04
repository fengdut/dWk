function f=xi_r(t,x)
global xs deltax
if x<=xs
    f=cos(t);
elseif x>=xs-deltax/2 && x<=xs+deltax/2
    f=(deltax-x+(xs-deltax/2)).*cos(t)./deltax;
else
    f=0.0;
end