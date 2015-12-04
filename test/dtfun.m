function f=dtfun(Lambda,t,x) % 1/(b*sqrt(1-Lambda/b))
f=1./(b(t,x).*sqrt(1-Lambda./b(t,x))); % best for simplist function of theta 
end
