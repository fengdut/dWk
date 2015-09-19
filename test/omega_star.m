function f=omega_star(epsilon,x)
global mn R0 a rhoh
f=mn*denp(x).*FE(epsilon)*R0*rhoh./...
    (a^2*2*x.*dFdE(epsilon).*density(x));
end
