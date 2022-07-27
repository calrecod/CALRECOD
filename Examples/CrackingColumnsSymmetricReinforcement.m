% CrackingColumnsSymmetricReinforcement
%----------------------------------------------------------------
% PURPOSE 
%    To design modify the momentum of inertia and effective concrete area
%    of a column cross-section symmetrically reinforced, considering a 
%    "Cracked" or "Non-Cracked" mechanisms of such cross-section due to 
%    load eccentricities
%
%    Note: widthEfficiencyCols function is used to obtain the itneraction
%          diagrams with respect to both cross-section axis directions, as
%          well as the neutral axis depth for each direction according to
%          the load combinations
%          
%          .diagramISR function is used to plot the interaction diagrams
%
%           CrackingColumnsSym function is the one used to compute the
%           modified inertia momentums of the cross-section for axis
%           directions, according the laod combinations (load eccentricity)
%           and neutral axis depth c
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

h=80;
b=60;
eccentricityXY=[120,115]; %cm
rec=[4,4];
fc=300; 
twidth=0.3;
height=300;
Pu=-41.3; %Ton

betac=1.05-fc/1400;
if betac<0.65
    betac=0.65;
elseif betac>0.85
    betac=0.85;
end
dimensionesColumna=[b,h];
Aco=b*h;
fdpc=0.85*fc;
Mux=abs(Pu*eccentricityXY(1)*0.01); Muy=abs(Pu*eccentricityXY(2)*0.01); % Ton-m
conditions=[1 Pu Mux Muy]; % Ton-m
E=2.1e6;

% Give reinforcement to the section 
t_widthISR=0.3;

% Plot the interaction diagram and load conditions for better appreciation

[Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(t_widthISR,...
        dimensionesColumna,rec,4200,30,conditions,fdpc,2e6,betac);


diagramISR(diagramaInteraccion,conditions);
% Compute, finally, the modified momentum of inertia for both axis 
% directions due to cracking mechanisms, as well as the transformed
% concrete area __________________________________________________
conditionCrack="Cracked";

Inertia0=[b*h^3/12,h*b^3/12]
tx=t_widthISR;
ty=t_widthISR;

% The neutral axis depths (cx,cy) must be determined based on the column's 
% interaction diagram in both directions and the load conditions: Take
% parameter [cxy] from the "widthEfficiencyCols" function;
[InertiaXYmodif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,...
          tx,eccentricityXY,ty,Pu,cxy,conditionCrack,E)