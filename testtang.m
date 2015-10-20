clear;
 
 
 n=500;
 
 
 delta=1e-10;
 omega=0.5;
 
 
 nd=250;
 
 delta=1e-10;
 yt=0;
 for ni=1:nd
     


 delta=delta*1.1;
 x(ni)=delta;
 yt(ni)=2*atan(1/(2*delta));
 
 

 end
 
 
 plot(x,yt);