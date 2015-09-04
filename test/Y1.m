function f=Y1(Lambda,t,tf,x,p)  % the integrand of Yp
%f=Gfun(Lambda,t,x).*exp(-i*p*tf)./...
%    (b(t,x)*sqrt(1-Lambda./b(t,x))); % tf: time as function of theta
f=exp(-i*p*tf);
end
