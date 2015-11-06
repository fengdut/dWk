function f=omega_star(denp,FE,dFdE,density,x)
global mn R0 a rhoh
f=mn*denp.*FE.*R0*rhoh./...
    (a^2*2*x.*dFdE.*density);
end
