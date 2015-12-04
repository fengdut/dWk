function f=denp(x) % derivative of density
global x0 deltax0
x1=(x-x0)./deltax0;
f=-2*x1.*exp(-x1.^2)./deltax0;
end
