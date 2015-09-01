function f=dWk(omega)
dx=0.1;
dE=0.5;
dL=0.1;
xa=0.0;
xb=1.0;
Ea=0.2;
Eb=4.0;
La=0;
nx=(xb-xa)/dx;
nxo2=nx/2;
nE=(Eb-Ea)/dE;
nEo2=nE/2;
tau=1.0;
sssum=0.0;
% using Simpson rule for three dimension integration
for j=1:nxo2
    
    x0=xa+dx*(j-1)*2.0
    x1=xa+dx*((j-1)*2.0+1);
    x2=xa+dx*((j-1)*2.0+2);
    ssum0=0.0;
    ssum1=0.0;
    ssum2=0.0;
    Lb=(1-eps1(x0)-0.1); % Lambda is a function of x
    nL=(Lb-La)/dL;
    nLo2=nL/2;
    for j1=1:nLo2
        
        L0=La+dL*(j1-1)*2.0;
        L1=La+dL*((j1-1)*2.0+1.0);
        L2=La+dL*((j1-1)*2.0+2.0);

        sum0=0.0;
        sum1=0.0;
        sum2=0.0;

        for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L0,epsilon0,x0,tau);
             ff1=WF(omega,L0,epsilon1,x0,tau);
             ff2=WF(omega,L0,epsilon2,x0,tau);
             sum0=sum0+dE/3*(ff0+4.0*ff1+ff2);
         end
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L1,epsilon0,x0,tau);
             ff1=WF(omega,L1,epsilon1,x0,tau);
             ff2=WF(omega,L1,epsilon2,x0,tau);
             sum1=sum1+dE/3*(ff0+4.0*ff1+ff2);
         end   
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L2,epsilon0,x0,tau);
             ff1=WF(omega,L2,epsilon1,x0,tau);
             ff2=WF(omega,L2,epsilon2,x0,tau);
             sum2=sum2+dE/3*(ff0+4.0*ff1+ff2);
         end
            ssum0=ssum0+dL/3*(sum0+4.0*sum1+sum2);
    end
    Lb=1-eps1(x1)-0.1;
    nL=(Lb-La)/dL;
    nLo2=nL/2;
    for j1=1:nLo2
        L0=La+dL*(j1-1)*2.0;
        L1=La+dL*((j1-1)*2.0+1.0);
        L2=La+dL*((j1-1)*2.0+2.0);
        sum0=0.0;
        sum1=0.0;
        sum2=0.0;
        for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L0,epsilon0,x1,tau);
             ff1=WF(omega,L0,epsilon1,x1,tau);
             ff2=WF(omega,L0,epsilon2,x1,tau);
             sum0=sum0+dE/3*(ff0+4.0*ff1+ff2);
         end
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L1,epsilon0,x1,tau);
             ff1=WF(omega,L1,epsilon1,x1,tau);
             ff2=WF(omega,L1,epsilon2,x1,tau);
             sum1=sum1+dE/3*(ff0+4.0*ff1+ff2);
         end
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L2,epsilon0,x1,tau);
             ff1=WF(omega,L2,epsilon1,x1,tau);
             ff2=WF(omega,L2,epsilon2,x1,tau);
             sum2=sum2+dE/3*(ff0+4.0*ff1+ff2);
         end
            ssum1=ssum1+dL/3*(sum0+4.0*sum1+sum2);
    end
    Lb=1-eps1(x2)-0.1;
    nL=(Lb-La)/dL;
    nLo2=nL/2;
    for j1=1:nLo2
        L0=La+dL*(j1-1)*2.0;
        L1=La+dL*((j1-1)*2.0+1.0);
        L2=La+dL*((j1-1)*2.0+2.0);
        sum0=0.0;
        sum1=0.0;
        sum2=0.0;
        for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L0,epsilon0,x2,tau);
             ff1=WF(omega,L0,epsilon1,x2,tau);
             ff2=WF(omega,L0,epsilon2,x2,tau);
             sum0=sum0+dE/3*(ff0+4.0*ff1+ff2);
         end
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L1,epsilon0,x2,tau);
             ff1=WF(omega,L1,epsilon1,x2,tau);
             ff2=WF(omega,L1,epsilon2,x2,tau);
             sum1=sum1+dE/3*(ff0+4.0*ff1+ff2);
         end
         for j2=1:nEo2
             epsilon0=Ea+dE*(j2-1)*2.0;
             epsilon1=Ea+dE*((j2-1)*2.0+1.0);
             epsilon2=Ea+dE*((j2-1)*2.0+2.0);
             ff0=WF(omega,L2,epsilon0,x2,tau);
             ff1=WF(omega,L2,epsilon1,x2,tau);
             ff2=WF(omega,L2,epsilon2,x2,tau);
             sum2=sum2+dE/3*(ff0+4.0*ff1+ff2);
         end
            ssum2=ssum2+dL/3*(sum0+4.0*sum1+sum2);
     end
        sssum=sssum+dx/3*(ssum0+4.0*ssum1+ssum2);
end
 f=sssum;
