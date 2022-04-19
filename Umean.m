% This function calculates the mean velocity and the friction factor
% This function receives:
% d,S,g,hc,CD,a,phi,alpha,beta,C,kappa,Pi as defined in Run.m
% This function returns:
% Depth-aeraged mean velocity (Umean); friction factor (f)

function [Umean,f] = Umean(d,S,g,hc,CD,a,phi,alpha,beta,C,kappa,Pi,UUD,y0)
%% Substitute the solved d to its relevant parameters
um=sqrt(g*S*(d-hc)); %shear velocity at the canopy top (m/s)
utot=sqrt(g*S*(d+hc*(phi-1))); %shear velocity at the channel bottom (m/s)
zeta=0.27*(utot/sqrt(2*g*S/(CD*a)))^(2.3)+1.15; %dimensionless inflection point velocity (zeta = ui/uUD)(-)
%% Calculation of the mean velocity of each layer
UUDmean = sqrt(2*g*S/(CD*a)); %uniform distribution layer
UMLmean = sqrt(2*g*S/(CD*a))*(zeta-1)*(1+(alpha*hc/d)*log(cosh(((beta/alpha)-d/(alpha*hc)))/cosh(beta/alpha))); %mixing layer
ULLmean = (1/d)*(um/kappa)*((d)-(beta*hc))*(log((d - beta*hc)/hc)+kappa*C-1); %log layer
UWFmean = um*Pi/kappa; %wake function layer

%% Calculation of the mean velocity
Umean = UUDmean+UMLmean+ULLmean+UWFmean; %summation of all the mean velocity components

%% Calculation of the friction factor
f=8*(utot/Umean)^2; 


end


