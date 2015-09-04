function f=pFpE(Lambda,epsilon) % the partial derivative of the distribution 
Lambda0=0.8;
dL=0.1;
f=-exp(-epsilon).*...
  exp(-(Lambda-Lambda0).^2/dL^2);
end
