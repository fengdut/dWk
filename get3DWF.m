

clear;

global R0 a rhoh 
R0=1.65;
a=0.4;
rhoh=0.1;
xa=0.0;
xb=1.0;

Ea=2.0;
Eb=4.0;

La=0.01;
Lb=1;

nx=22;
nL=13;
nE=16; 



dx=(xb-xa)/(nx-1);
dL=(Lb-La)/(nL-1);
dE=(Eb-Ea)/(nE-1);

tau=1.0;


omega=0.1+0.01i;

WF3D(1:nx,1:nL,1:nE)=0;

for  jx=2:nx
    x=xa+(jx-1)*dx;
    jx
    for jL=1:nL
        L = La+(jL-1)*dL;      
        for jE=1:nE
          E=Ea + (jE-1)*dE;
          WF3D(jx,jL,jE)=WF(omega,L,E,x,tau);
        end
     end
  
end


sum2d(1:nx,1:nL)=0;
% https://en.wikipedia.org/wiki/Simpson%27s_rule  Simpson's 3/8 rule (for n intervals)
for  jx=1:nx
    for jL=1:nL
        
         for jE=1:nE
             sum2d(jx,jL)=sum2d(jx,jL)+ WF3D(jx,jL,jE);
         end
         for j1=2:nE-1
             sum2d(jx,jL)=sum2d(jx,jL)+ 2*WF3D(jx,jL,jE);
         end
         for j1=4:3:nE-1
             sum2d(jx,jL)=sum2d(jx,jL)- WF3D(jx,jL,jE);
         end
    end
end

sum2d=sum2d*3*dE/8;


sum1d(1:nx)=0;
for  jx=1:nx
    
    for jL=1:nL
    sum1d(jx) =sum1d(jx) +sum2d(jx,jL);
    end
    for jL=2:nL-1
    sum1d(jx) =sum1d(jx) +2*sum2d(jx,jL);
    end
    for jL=4:3:nL-1
    sum1d(jx) =sum1d(jx)-sum2d(jx,jL);
    end
end
sum1d=sum1d*3*dL/8;

sumall=0;
for  jx=1:nx
    sumall = sumall+ sum1d(jx);
end

for  jx=2:nx-1
    sumall = sumall+ 2*sum1d(jx);
end

for  jx=4:3:nx-1
    sumall = sumall- sum1d(jx);
end
sumall=sumall*3*dx/8



%ntt=nx*nL*nE;

%WF1D=reshape(WF3D,ntt,1);
%sum(sum(sum(WF3D)))
% 
% parfor j=1:ntt
%     
%     jx=ceil(j/(nL*nE))
%     rntt = j -(jx-1)*(nL*nE);
%     jL=ceil(rntt/nE)
%     jE=rntt -(jL-1)*nE
%     
%     x=xa+(jx-1)*dx;
%     L = La+(jL-1)*dL;
%     E=Ea + (jE-1)*dE;
%     tx=WF(omega,L,E,x,tau)
%     WF1D(j)=0;
%     
% end
