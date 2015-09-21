function f=xi_r(t,x)
global xs deltax xi0 nt
f=zeros(1,nt+1);
for j=1:nt+1
if x(j)<=xs
    f(j)=xi0*cos(t(j)); % displacement varies with theta
elseif x(j)>=xs-deltax/2 && x(j)<=xs+deltax/2
    f(j)=xi0*(deltax-x(j)+(xs-deltax/2)).*cos(t(j))./deltax;
else
    f(j)=0.0;
end
end
end
