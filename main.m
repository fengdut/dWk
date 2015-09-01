clear;
close all;
global R0 a rhoh 
R0=1.65;
a=0.4;
rhoh=0.1;
dwr=0.01;
wra=0.01;
wrb=4.01;
nr=(wrb-wra)/dwr;
dwi=0.01;
wia=0.01;
wib=0.51;
ni=(wib-wia)/dwi;
dWF=0.5; % MHD contribution to delta W
w=[];
dWk(0.1+0.01*i)
%f=@(x)i*x+dWF+dWk(x); % fishbone dispersion relation, iomega/omega_A+dWF+dWk(omega)=0
%for j=1:ni
%    wi=wia+(j-1)*dwi;
%    for j1=1:nr
%        wr=wra+(j1-1)*dwr;
%        options=optimset('Display','off');
%        x=fsolve(f,wr+i*wi,options);
%        w=[w,x];
%    end
%end
