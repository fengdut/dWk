function [F]=simpintegral_3D(F_3D,nx,dx,nL,dL,nE,dE)

F_2D(1:nL,1:nE)=0;


    for iL=1:nL
        for iE=1:nE
        F_2D(iL,iE)= simpintegral(F_3D(:,iL,iE),nx,dx);
        end
    end


F_1D(1:nL)=0;

for iE=1:nE
    F_1D(iE)= simpintegral(F_2D(:,iE),nL,dL);
end

F=simpintegral(F_1D,nE,dE);
end