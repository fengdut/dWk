function f=xi_t(t,x)
global xs deltax
if x<=xs
    f=-x*sin(t);
elseif x>=xs-deltax/2 && x<=xs+deltax/2
    f=-(deltax-x+(xs-deltax/2)).*x.*sin(t)./...
        deltax-x.^2.*sin(t)/deltax;
else
    f=0.0;
end