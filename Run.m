%% Solver for flow resistance on vegetated channels 
%Version 1.0 (April 2022)
%Developed by: Matthias Kramer and Hanwen Cui
%Works with Matlab 2021a
%Contact:
%m.kramer@adfa.edu.au

%When using this code, please cite/refer to the following reference:
%--------------------------------------------------------------------------%
%H. Cui, S. Felder, and M. Kramer (2022)
%Predicting flow resistance in vegetated open-channel flows 
%Water Resources Research
%--------------------------------------------------------------------------%

clear all
close all

%% Input parameters
q=0.05; %define the flow rate that needs to be calculated
X = ['q (input) = ', num2str(q),' (m^2/s)']; disp (X) 
syms d %water depth (m)

% Channel characteristics
S=sind(10.8); %chute slope (-)
g=9.81; %(m/s^2)

% Vegetation properties
hc=0.01; %deflected canopy height (m)
CD=1.2; %drag coefficient (-)
a=10; %frontal area per volume (1/m)
phi=1; %porosity (-)

%% Classify vegetation type
Den=CD*a*hc; %hydrodynamic vegetation density (-)

% Dense canopy
if Den > 0.5 % boundary for dense canopy
    beta = 0.85; %dimensionless elevation of the inflection point (beta=yi/hc) (-)
    C=4.0; %integration constant log-law (-)
    Pi=0.3; %wake function term (-)
else %Transitional canopy
    beta = 0.77; %dimensionless elevation of the inflection point (beta=yi/hc) (-)
    C=2.2; %integration constant log-law (-)
    Pi=0.15; %wake function term (-)
end

%% Model parameters
alpha=0.1+0.01/Den; %dimensionless mixing layer length scale (alpha=Le/hc) (-)
kappa=0.41; %von Karman constant (-)
UUD=sqrt(2*g*S/(CD*a)); %inside-canopy uniform velocity;
y0=hc*exp(-kappa*C); %hydraulic roughness(m)
um=sqrt(g*S*(d-hc)); %shear velocity at the canopy top (m/s)
utot=sqrt(g*S*(d+hc*(phi-1))); %shear velocity at the channel bottom (m/s)
zeta=0.27*(utot/sqrt(2*g*S/(CD*a)))^(2.3)+1.15; %dimensionless inflection point velocity (zeta = ui/uUD)(-)

%% Numerical Solution for d
um=sqrt(g*S*(d-hc)); %shear velocity at the canopy top(m/s)
eqnLeft=q/d; %continuity equation
eqnRight=sqrt(2*g*S/(CD*a))+sqrt(2*g*S/(CD*a))*(zeta-1)*(1+(alpha*hc/d)*log(cosh(((beta/alpha)-d/(alpha*hc)))/cosh(beta/alpha)))+(1/d)*(um/kappa)*((d)-(beta*hc))*(log((d - beta*hc)/hc)+kappa*C-1)+um*Pi/kappa; %Eq. (2) in ......

% Numerical solver for d
Sol = vpasolve(eqnLeft == eqnRight, d);
dSol= double(Sol);

% Calculation of the mean velocity and friction factor
[UmSol,fSol]=Umean(dSol,S,g,hc,CD,a,phi,alpha,beta,C,kappa,Pi,UUD,y0);
dSolhc=dSol/hc; %submergence of the interested grass-lined flow
FrSol=UmSol/sqrt(dSol*g); % Froude number of the flow

X = ['d (solution) = ', num2str(dSol),' (m)']; disp (X) 
X = ['u (solution) = ', num2str(UmSol),' (m/s)']; disp (X) 
X = ['f (solution) = ', num2str(fSol), ' (-)']; disp (X) 
X = ['Fr (solution) = ', num2str(FrSol), ' (-)']; disp (X) 