function f=density(x) % density profile
global x0 deltax0
x1=(x-x0)./deltax0;
f=exp(-x1.^2);
end
