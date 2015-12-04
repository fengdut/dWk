clear
close all
global R0 a rhoh mn nn e xs deltax
global nL nx nt ne xa xb La Lb ta tb ea eb np
global dt
global p rhod xi0
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

omega=1.5+0.1e-02i;

mn=1; % poloidal mode number
nn=1; % toroidal mode number

nL=50; % Lambda grid number 
nt=50; % theta grid number for poloidal angle
nx=50; % x grid number for position
ne=50; % epsilon grid number for energy


xa=1e-06; % left boundary of x
xb=1.0;   % right boundary of x
La=1e-06; % left boundary of Lambda
Lb=1-0.3; % right boundary of Lambda
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
    Lb=1-eps1(x(j));           
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

tfun3D=zeros(nL+1,nx+1,nt+1); % collect tfun(Lambda,x,theta) data 
                              % store as 3D array
                              % according to above 
                              % defined coordinate grids


p=-1:1;                       % p=-1,0,1
np=length(p);                 

Yp4D=zeros(np,ne+1,nL+1,nx+1); % collect Yp(p,epsilon,Lambda,x) data 
                          %  store as 3D array
                          %  according to above 
                          %  defined coordinate grids
                          
WF3D=zeros(ne+1,nL+1,nx+1); % collect WF(epsilon,Lambda,x) data 
                            % store as 3D array
                            % according to above 
                            % defined coordinate grids

tic
for j=1:nL                % j=nL+1 ellipitic function become infinite.
    for j1=1:nx+1         % elliminte this point
        for j2=1:nt+1
        tfun3D(j,j1,j2)=tfun(Lambda(j,j1),t(j2),x(j1),tau); % the data of 
        end                             % tfun store as 3D array, 
                                        % which also accelerates the 
                                        % speed without calculating 
                                        % tfun every step
        
    end
end
toc
tic
for j2=1:np
 for j3=1:ne+1
  for j=1:nL
    for j1=1:nx+1
        Yp4D(j2,j3,j,j1)=Yp(Lambda(j,j1),ee(j3),tfun3D(j,j1,:),x(j1),p(j2),tau); % the data 
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
for j2=1:ne+1
  for j=1:nL
    for j1=1:nx+1
        WF3D(j2,j,j1)=WF(omega,Lambda(j,j1),ee(j2),x(j1),Yp4D(:,j2,j,j1)); % the data 
    end                                      % of WF store as 3D 
                                             % array, which also 
                                             % accelerates the speed 
                                             % without calculating WF 
                                             % every step
        
    
  end
end
toc
% calculating delta W_k 3D integral
tic
WFe=zeros(1,ne+1);
I1=zeros(1,nL);
I2=zeros(1,nx+1);

for j1=1:nx+1
    for j=1:nL
        WFe(:)=WF3D(:,j,j1);
        I1(j)=simp(ne+1,de,WFe); % integral of WF with 
                                 % respect to epsilon at 
                                 % fixed x and for nL Lambda's
    end
    I2(j1)=simp(nL,dL,I1); % integral of I1 with 
                           % respect to Lambda for
                           % nx+1 x's
end
I3=simp(nx+1,dx,I2); % intergral of I2 with 
                     % respect to x
toc                     
                     
                     