function f=FE(epsilon) % distribution of energy compenont
global epsilonc epsilon0 deltae
epsilon1=(epsilon-epsilon0)/deltae;
f=1./(epsilon.^1.5+epsilonc.^1.5).*...
    erf(epsilon1);
end
