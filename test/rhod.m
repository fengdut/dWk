function f=rhod(Lambda,epsilon,t,x) % drift orbit width for passing particles
global rhoh
f=qprofile(x).*rhoh/2.*sqrt(epsilon./(1-Lambda./b(t,x))).*...
    (Lambda./b(t,x)+2*(1-Lambda./b(t,x)));
end