function f=xi_r(t,x)
global xs deltax xi0 nt
for j=1:nt+1
if x(j)<=xs
    f=xi0*cos(t);
elseif x(j)>=xs-deltax/2 && x(j)<=xs+deltax/2
    f=xi0*(deltax-x(j)+(xs-deltax/2)).*cos(t)./deltax;
else
    f=0.0;
end
end
end
