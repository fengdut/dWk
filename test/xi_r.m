function f=xi_r(t,x)
global xs deltax xi0
if x<=xs
    f=xi0*cos(t); % displacement varies with theta
elseif x>=xs-deltax/2 && x<=xs+deltax/2
    f=xi0*(deltax-x+(xs-deltax/2)).*cos(t)./deltax;
else
    f=0.0;
end
end
