 clear;
 
 
 n=50000;
 
 
 delta=1e-10;
 omega=0.5;
 
 
 nd=250;
 
 for ni=1:nd
     
 omegad=linspace(0.0,1.0,n);
 
 
 %omegad=omegad*
 delta=delta*1.1;
 
 
 f=delta./((omega-omegad).^2+delta^2);
 
 deltaa(ni)=delta;
 sumf(ni)=simpintegral(f,n,omegad(2)-omegad(1));
 
 end
 
 
 plot(deltaa,(sumf),'o--');
 
 %hold all;
 %plot(deltaa,1./deltaa,'+');
 
 
 
 
 
 
 
 
%  x=linspace(0,pi,n);
%  ycos=sin(x);
%  
%  %plot(x,ycos);
%  x(2)-x(1)
%  sy=simpintegral(ycos,n,x(2)-x(1))

 
 
