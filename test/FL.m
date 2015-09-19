function f=FL(Lambda) % distribution of Lambda component
deltaL=0.8;
L0=0.1;
f=exp(((Lambda-L0)./deltaL^0.2).^2);
end
