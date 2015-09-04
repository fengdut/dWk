clear;

tn=100;
dx=(2*pi)/(tn-1);
x=0:dx:2*pi;

y=sin(x);


sumy=simpintegral(y,tn,dx)
quad(@(x)sin(x),0,2*pi)

