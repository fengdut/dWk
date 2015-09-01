function f=kappa_r(t,x)
global e
e=0.33;
f=-e*cos(t)-e.^2./qprofile(x);
end
