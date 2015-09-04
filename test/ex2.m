clear
close all
nL=20; % Lambda grid number
nt=500; % theta grid number
nx=10; % x grid number
ne=20; % epsilon grid number

xa=1e-06; % left boundary of x
xb=1.0;   % right boundary of x
La=1e-06; % left boundary of Lambda
Lb=1-0.3; % right boundary of Lambda
ta=1e-06; % left boundary of theta
tb=2*pi;  % right boundary of theta
ea=1e-06; % left boundary of epsilon
eb=4.0;   % right boundary of epsilon

dx=(xb-xa)/(nx);
x=xa:dx:xb; % construct x array priorly, 
            % without calculating xi every step
            
Lambda=zeros(nL+1,length(x)); 
for j=1:length(x)
    Lb=1-eps1(x(j));
    dL=(Lb-La)/(nL);
    Lambda(:,j)=La:dL:Lb;
end