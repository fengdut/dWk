function f=tfun(Lambda,epsilon,t,x,tau)

tn=size(t,2);
f(1:tn)=0;

for i=2:tn
    
f(i)=quad(@(x)dtfun(Lambda,epsilon,t(i),x,tau),0,t(i));

end
end