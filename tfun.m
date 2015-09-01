function f=tfun(Lambda,epsilon,t,x,tau)

tn=size(t,2);
f(1:tn)=0;

for i=1:tn
    
%f(i)=quad(@(t)dtfun(Lambda,epsilon,t,x,tau),0,t(i));

%f(i)=quad(@(t)qprofile(x)./sqrt(2*epsilon)./b(t,x)./sqrt(1-Lambda./b(t,x))*tau,0,t(i));

f(i)=quad(@(t)qprofile(x)./(sqrt(2*epsilon).*b(t,x).*sqrt(1-Lambda./b(t,x)).*tau),0,t(i));

%Fun=@(t)dtfun(Lambda,epsilon,t,x,tau);

%f(i) = integeral(Fun,0,t(i));
end
end