function f=Y1(Gfun,tfun,p,Lambda,b)  % the integrand of Yp
f=Gfun.*exp(-i*p.*tfun)./...
    (b.*sqrt(1-Lambda./b)); % tf: time as function of theta
%f=exp(-i*p*tf);                      % relates to the corresponding theta
end
