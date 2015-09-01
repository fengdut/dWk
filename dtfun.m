function f=dtfun(Lambda,epsilon,t,x,tau)
 f=qprofile(x)./sqrt(2*epsilon)./...
     b(t,x)./sqrt(1-Lambda./b(t,x))*tau;


end