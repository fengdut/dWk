function f=xi_t(t,x)
global xs deltax xi0 
if x<=xs
    f=-xi0*x.*sin(t); % displacement varies with theta
elseif x>=xs-deltax/2 && x<=xs+deltax/2
    f=xi0*(-(deltax-x+(xs-deltax/2)).*x.*sin(t)./...
        deltax+x.^2.*sin(t)/deltax);
else
    f=0.0;
end
end
