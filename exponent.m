function f=exponent(Lambda,epsilon,t,x,tau,p)
f=exp(-i*p*omega_b(Lambda,epsilon,x)*...
       tfun(Lambda,epsilon,t,x,tau));
end