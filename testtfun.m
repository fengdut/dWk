

global R0 a rhoh xa xb La Lb Ea Eb dx dL dE
global nx nL nE 

global tfun3d
nt=10;


ta=0;
tb=2*pi;
dt=(tb-ta)/(nt-1);

tfun3d(1:nx,1:nL,1:nt)=0;

for ix=1:nx
    x=xa+(ix-1)*dx;
    for il=1:nL
        L=La+dL*(il-1);
        for it=2:nt
            t=ta+(it-1)*dt;
            if(L>=1-eps1(x))
                
            else
             tfun3d(ix,il,it) = tfun3d(ix,il,it)+ dt*(1+x^2)/((1-eps1(x)*cos(t))*sqrt(1-L./(1-eps1(x)*cos(t))));
            end
        end
    end
end

