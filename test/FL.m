function f=FL(Lambda) % distribution of Lambda component
global Lambda0 deltaL
Lambda1=(Lambda-Lambda0)./deltaL.^0.2;
f=exp(Lambda1.^2);
end
