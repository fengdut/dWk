function f=rhod_x(qprofile,b,Lambda,epsilon,t,x) % calculate x+rhod*cos(t)
global rhoh
f=x+qprofile.*rhoh/2.*sqrt(epsilon./(1-Lambda./b)).*...
    (Lambda./b+2*(1-Lambda./b)).*cos(t);
end
