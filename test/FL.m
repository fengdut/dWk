function f=FL(Lambda) % distribution of Lambda component
global Lambda0 deltaL
Lambda1=(Lambda-Lambda0)./deltaL;
f=exp(-Lambda1.^2);
end
