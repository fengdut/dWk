function f=xi_t(t,x)
global xs deltax xi0 nt
for j=1:nt+1
if x(j)<=xs
    f=-xi0*x(j).*sin(t);
elseif x(j)>=xs-deltax/2 && x(j)<=xs+deltax/2
    f=xi0*(-(deltax-x(j)+(xs-deltax/2)).*x(j).*sin(t)./...
        deltax+x(j).^2.*sin(t)/deltax);
else
    f=0.0;
end
end
end
