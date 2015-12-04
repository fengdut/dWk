clear
close all
global R0 a rhoh mn nn e xs deltax
global nL nx nt ne xa xb La Lb ta tb ea eb np
global dt tau
global p xi0
global epsilonc epsilon0 deltae
global x0 deltax0
global Lambda0 deltaL

R0=1.65; % major radius unit m
a=0.40;  % minor radius
e=a/R0;
rhoh=0.08; % vh/Omega=sqrt(2Th/M)/Omega, Omega=Be/M
xs=0.5; % rs/a position of rational surface q=1
deltax=0.06; % the artificial width along q=1 surface
xi0=1.0; % ratio of displacement to minor radius, xi0/a

% omega 
omega_min=1.5; 
omega_max=1.5; 
domega=0.1;
omega=omega_min:domega:omega_max;
no=length(omega);
omega=omega+0.1e-02i; 

mn=1; % poloidal mode number
nn=1; % toroidal mode number

nL=22; % Lambda grid number 
nt=52; % theta grid number for poloidal angle
nx=52; % x grid number for position
ne=22; % epsilon grid number for energy


xa=1e-06; % left boundary of x
xb=1.0;   % right boundary of x
La=1e-06; % left boundary of Lambda
%Lb=1-0.3; % right boundary of Lambda
ta=1e-06; % left boundary of theta
tb=2*pi;  % right boundary of theta
ea=1e-01; % left boundary of epsilon
eb=5.0;   % right boundary of epsilon

% distribution function
epsilonc=0.8*eb; % critical energy of slow down distribution
epsilon0=eb; % 
deltae=1.0; % energy width
x0=0.0; % injection position of NBI
deltax0=0.3; % density width
Lambda0=0.05; % injection Lambda of NBI
deltaL=0.02; % Lambda density width


% coordinate grids

dx=(xb-xa)/(nx);
x=xa:dx:xb; % construct x array priorly, 
            % without calculating xi every step
            
Lambda=zeros(nL+1,length(x));  % construct Lambda 2D array due to
for j=1:length(x)              % Lambda maximum as a function of x  
    %Lb=1-eps1(x(j)); 
    Lb=0.1;
    dL=(Lb-La)/(nL);
    Lambda(:,j)=La:dL:Lb;
end
            

dt=(tb-ta)/nt;
t=ta:dt:tb;  %left boundary of theta  should be finite, 
             % otherwise, cannot form the grids of t for caculating theta
             % array t construct priorly which will 
             %accelerate the speed without caculating tj every step.
             
de=(eb-ea)/ne;
ee=ea:de:eb;

tau=1.0; % +1 or -1 for v_parallel direction 

p=-1:1;                       % p=-1,0,1
np=length(p);

qprofile1D=zeros(1,nx+1);
Jac1D=zeros(1,nx+1);
density1D=zeros(1,nx+1);
dFdE1D=zeros(1,ne+1);
FE1D=zeros(1,ne+1);
denp1D=zeros(1,nx+1);
eps1_1D=zeros(1,nx+1);
FL2D=zeros(nL+1,nx+1);
qprofile1D(:)=qprofile(x);
Jac1D(:)=Jac(x);
density1D(:)=density(x);
denp1D(:)=denp(x);
dFdE1D(:)=dFdE(ee);
FE1D(:)=FE(ee);
eps1_1D(:)=eps1(x);



b2D=zeros(nt+1,nx+1);
kappa2D=zeros(nL+1,nx+1);
K2D=zeros(nL+1,nx+1);
gtt2D=zeros(nt+1,nx+1);
grt2D=zeros(nt+1,nx+1);
grr2D=zeros(nt+1,nx+1);
kappa_t2D=zeros(nt+1,nx+1);
kappa_r2D=zeros(nt+1,nx+1);
omega_star2D=zeros(ne+1,nx+1);
omega_b3D=zeros(nL+1,ne+1,nx+1);
omega_phi3D=zeros(nL+1,ne+1,nx+1);
omega_minus_omega_star3D=zeros(no,ne+1,nx+1);
R_term5D=zeros(np,no,nL+1,ne+1,nx+1);
fR_term1D=zeros(1,np);
rhod_x4D=zeros(nL+1,ne+1,nt+1,nx+1);
Gfun4D=zeros(nL+1,ne+1,nt+1,nx+1);
xi_t4D=zeros(nL+1,ne+1,nt+1,nx+1);
xi_r4D=zeros(nL+1,ne+1,nt+1,nx+1);
fxi_r1D=zeros(1,nx+1);
fxi_t1D=zeros(1,nx+1);
Y1_5D=zeros(np,nL+1,ne+1,nt+1,nx+1);
fY1_1D=zeros(1,nt+1);


tfun3D=zeros(nL+1,nt+1,nx+1); % collect tfun(Lambda,x,theta) data 
                              % store as 3D array
                              % according to above 
                              % defined coordinate grids


                 

Yp4D=zeros(np,ne+1,nL+1,nx+1); % collect Yp(p,epsilon,Lambda,x) data 
fYp1D=zeros(1,np);                          %  store as 3D array
                          %  according to above 
                          %  defined coordinate grids
                          
WF4D=zeros(no,ne+1,nL+1,nx+1); % collect WF(epsilon,Lambda,x) data 
                            % store as 3D array
                            % according to above 
                            % defined coordinate grids

tic

for j=1:nL

    for j1=1:nx+1
        kappa2D(j,j1)=kappa(Lambda(j,j1),eps1_1D(j1));
        FL2D(j,j1)=FL(Lambda(j,j1));
        K2D(j,j1)=ellipke(1./kappa2D(j,j1));
     end
end
toc
tic
for j=1:ne+1

    for j1=1:nx+1
        omega_star2D(j,j1)=omega_star(denp1D(j1),FE1D(j),dFdE1D(j),...
            density1D(j1),x(j1));
     end
end
toc
tic
for j2=1:no
  for j1=1:ne+1

     for j=1:nx+1
        omega_minus_omega_star3D(j2,j1,j)=...
            real(omega(j2))-omega_star2D(j1,j);
     end
  end
end
toc
tic
for j=1:nL
    for j1=1:ne+1

       for j2=1:nx+1
        omega_b3D(j,j1,j2)=omega_b(kappa2D(j,j2),qprofile1D(j2),...
            Lambda(j,j2),ee(j1),eps1_1D(j2),K2D(j,j2));
        omega_phi3D(j,j1,j2)=qprofile1D(j2).*omega_b3D(j,j1,j2);
       end
    end
end
toc
tic
for j4=1:np
 for j3=1:no
  for j=1:nL
    for j1=1:ne+1

       for j2=1:nx+1
        R_term5D(j4,j3,j,j1,j2)=p(j4).*omega_b3D(j,j1,j2)+...
        nn.*omega_phi3D(j,j1,j2)-omega(j3);
       end
    end
  end
 end
end
toc
tic
for j1=1:nt+1
    for j=1:nx+1
        gtt2D(j1,j)=gtt(t(j1),x(j));
        grt2D(j1,j)=grt(t(j1),x(j));
        grr2D(j1,j)=grr(t(j1),x(j));
        kappa_r2D(j1,j)=kappa_r(t(j1),x(j));
        kappa_t2D(j1,j)=kappa_t(t(j1),x(j));
        b2D(j1,j)=b(t(j1),x(j));
    end
end
toc

tic
for j3=1:nL
    for j2=1:ne+1
        for j1=1:nt+1
            %for j=1:nx+1
                rhod_x4D(j3,j2,j1,:)=rhod_x(qprofile1D,...
                    b2D(j1,:),Lambda(j3,:),ee(j2),t(j1),x);
            %end
        end
    end
end
toc
tic
for j3=1:nL
    for j2=1:ne+1
        for j1=1:nt+1
            for j=1:nx+1
                xi_r4D(j3,j2,j1,j)=xi_r(t(j1),rhod_x4D(j3,j2,j1,j));
                xi_t4D(j3,j2,j1,j)=xi_t(t(j1),rhod_x4D(j3,j2,j1,j));
            end
        end
    end
end
toc
tic
for j3=1:nL
    for j2=1:ne+1
        for j1=1:nt+1
            %for j=1:nx+1
                fxi_r1D(:)=xi_r4D(j3,j2,j1,:);
                fxi_t1D(:)=xi_t4D(j3,j2,j1,:);
                Gfun4D(j3,j2,j1,:)=Gfun(gtt2D(j1,:),grt2D(j1,:),grr2D(j1,:),...
                    kappa_t2D(j1,:),kappa_r2D(j1,:),...
                    Lambda(j3,:),b2D(j1,:),...
                    fxi_t1D,...
                    fxi_r1D);
            %end
        end
    end
end
toc


tic
for j=1:nL                % j=nL+1 ellipitic function become infinite.
    for j1=1:nt+1         % elliminte this point
        for j2=1:nx+1
        tfun3D(j,j1,j2)=tfun(kappa2D(j,j2),Lambda(j,j2),t(j1),x(j2),...
            eps1_1D(j2),K2D(j,j2)); % the data of 
        end                             % tfun store as 3D array, 
                                        % which also accelerates the 
                                        % speed without calculating 
                                        % tfun every step
        
    end
end
toc
tic
for j4=1:np
  for j3=1:nL
    for j2=1:ne+1
        for j1=1:nt+1
            for j=1:nx+1
                Y1_5D(j4,j3,j2,j1,j)=Y1(Gfun4D(j3,j2,j1,j),...
                    tfun3D(j3,j1,j),p(j4),Lambda(j3,j),b2D(j1,j));
            end
        end
    end
  end
end
toc
tic
for j2=1:np
 for j3=1:ne+1
  for j=1:nL
    for j1=1:nx+1
        fY1_1D(:)=Y1_5D(j2,j,j3,:,j1);
        Yp4D(j2,j3,j,j1)=Yp(fY1_1D,...
            kappa2D(j,j1),Lambda(j,j1),eps1_1D(j1),K2D(j,j1)); % the data 
    end                                      % of Yp store as 3D 
                                              % array, which also acc-
                                              % -elerates the speed 
                                            % without calculating Yp 
                                            % every step
        
    
  end
 end
end
toc
tic
for j3=1:no
 for j2=1:ne+1
  for j=1:nL
    for j1=1:nx+1
        fYp1D(:)=Yp4D(:,j2,j,j1);
        fR_term1D(:)=R_term5D(j3,:,j,j2,j1);
        WF4D(j3,j2,j,j1)=WF(fYp1D,...
            fR_term1D,...
            Jac1D(j1),qprofile1D(j1),...
            ee(j2),density1D(j1),dFdE1D(j2),...
            FL2D(j,j1),omega_b3D(j,j2,j1),...
            omega_minus_omega_star3D(j3,j2,j1)); % the data 
    end                                      % of WF store as 3D 
                                             % array, which also 
                                             % accelerates the speed 
                                             % without calculating WF 
                                             % every step
        
    
  end
 end
end
toc
% calculating delta W_k 3D integral
tic
I3=zeros(1,no);
for j2=1:no
   WFe=zeros(1,ne+1);
   I1=zeros(1,nL);
   I2=zeros(1,nx+1);

 for j1=1:nx+1
    for j=1:nL
        WFe(:)=WF4D(j2,:,j,j1);
        I1(j)=simp(ne+1,de,WFe); % integral of WF with 
                                 % respect to epsilon at 
                                 % fixed x and for nL Lambda's
    end
    I2(j1)=simp(nL,dL,I1); % integral of I1 with 
                           % respect to Lambda for
                           % nx+1 x's
 end
    I3(j2)=simp(nx+1,dx,I2); % intergral of I2 with 
                     % respect to x
end
toc

                     
                     
