function f=omega_star(denp,FE,dFdE,density,x)
global mn e rhoh
f=mn*denp.*FE.*rhoh./...
    (e*2*x.*dFdE.*density);
end
