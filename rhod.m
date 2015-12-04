function f=rhod(qprofile,b,Lambda,epsilon) % drift orbit width for passing particles
global rhoh
f=qprofile.*rhoh/2.*sqrt(epsilon./(1-Lambda./b)).*...
    (Lambda./b+2*(1-Lambda./b));
end
