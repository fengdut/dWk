function f=kappa_r(t,x)
global e
f=-e*cos(t)./(1+e*x.*cos(t));
end
