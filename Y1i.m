%function f=Y1(Lambda,epsilon,t,x,tau,p)  % the integrand of Yp

% f=Gfun(Lambda,t,x).*qprofile(x)/...
%        sqrt(2*epsilon)./tau./b(t,x)./...
%         sqrt(1-Lambda./b(t,x)) .*...
%         exponent(Lambda,epsilon,t,x,tau,p);


%bt=b(t,x);
e=0.23;
bt=1-e*x*cos(t);
global tomega_b
global Y13D
global R0 a rhoh xa xb La Lb Ea Eb dx dL dE
global nx nL nE 


Y13D(1:nx,1:nL,1:nE)=0;



% f=Gfun(Lambda,t,x).*qprofile(x)./...
%        sqrt(2*epsilon)./tau./bt./...
%         sqrt(1-Lambda./bt).*exp(-i*p*tomega_b*tfun(Lambda,epsilon,t,x,tau));
    
    


 f=Gfun(Lambda,t,x).*qprofile(x)./...
    sqrt(2*epsilon)./tau./bt./...
    sqrt(1-Lambda./bt).*exp(-i*p*tomega_b*tfun3d(ix,iL,it))/(sqrt(2*epsilon));


for  jx=1:nx
    x=xa+(jx-1)*dx;

    for jL=1:nL
        L = La+(jL-1)*dL;      
        for jE=1:nE
          E=Ea + (jE-1)*dE;
          Y13D(jx,jL,jE)=Gfun(Lambda,t,x).*qprofile(x)./...
                          sqrt(2*epsilon)./tau./bt./...
                          sqrt(1-Lambda./bt).*exp(-i*p*tomega_b*tfun3d(ix,iL,it))/(sqrt(2*epsilon));

        end
     end
  
end

%end
