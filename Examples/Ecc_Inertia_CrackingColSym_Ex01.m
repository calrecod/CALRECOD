% Ecc_Inertia_CrackingColSym_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the modification of the inertia of a rectangular column
%    cross-section subject to biaxial bending-compression forces as the
%    load eccentricities increase
%
%    Note: widthEfficiencyCols function is used to obtain the interaction
%          diagrams with respect to both cross-section axis directions, as
%          well as the neutral axis depth for each direction according to
%          the load combinations
%          
%           CrackingColumnsSym function is the one used to compute the
%           modified inertia momentums of the cross-section for axis
%           directions, according the laod combinations (load eccentricity)
%           and neutral axis depth c
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-05-02
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
h=80; % cross-section's height dimension (cm)
b=60; % cross-section's width dimension (cm)
rec=[4,4]; % concrete cover for each cross-section direction (cm)
height=300; % column's length (cm)
dimensionsColumn=[b,h];

%% Materials
fc=300; % concrete's compressive strength (Kg/cm2)
fdpc=0.85*fc;
betac=1.05-fc/1400;
if betac<0.65
    betac=0.65;
elseif betac>0.85
    betac=0.85;
end
fy=4200; % yield stress of the reinforcing steel (Kg/cm2)
Es=2.1e6; % Modulus of Elasticity of the reinforcing steel (Kg/cm2)

%% Reinforcement data
% ISR's width to propose a reinforcement area
twidthISR=0.45; % cm

%% Additional computations
Aisr=(2*(b-rec(1))+2*(h-rec(2)))*twidthISR % proposed reinforcement area
Amin=0.01*b*h % Min reinforcement area by code
Amax=0.04*b*h % Max reinforcement area by code

Inertia0=[b*h^3/12,h*b^3/12]; % gross section's inertia

%% Loads
Pu=-105.3e3; % Compression load (Kg)

%% Additional parameters
conditionCrack="Cracked";

%% Main loop
IxyredVec=[];
eccxyVec=[];
maxEc=100; % Max load eccentricity for both section's axis
ecxy=0; % Initial load eccentricity
while ecxy<=maxEc
    ecxy=ecxy+0.1;
    eccentricityXY=[ecxy,ecxy]; % Load eccentricities for each 
                              % section's axis (cm)
    eccxyVec=[eccxyVec; eccentricityXY]; % collecting eccentricities
                                         % for plot
    Mux=abs(Pu*eccentricityXY(1)); % bending loads Mu = P * e
    Muy=abs(Pu*eccentricityXY(2)); % Kg-cm

    conditions=[1 Pu Mux Muy]; % Kg-cm

    %% Structural efficiency
    [Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols...
       (twidthISR,dimensionsColumn,rec,fy,60,conditions,fdpc,Es,betac);

    %% Modified momentum of inertia for both axis directions
    % cracking mechanisms are considered through the 
    % "transform section" method

    % Take parameter [cxy] from the "widthEfficiencyCols" function
    [InertiaXYmodif,Atransfxy,elimxy]=CrackingColumnsSym(h,b,fdpc,rec,...
              twidthISR,eccentricityXY,twidthISR,Pu,cxy,conditionCrack,Es);

    %% Computation of the inertia reduction degree
    redIxy=InertiaXYmodif./Inertia0;
    IxyredVec=[IxyredVec; redIxy];
end

%% Plotting results
figure(1)
plot(eccxyVec(:,1), IxyredVec(:,1)*100,'k -','LineWidth',1.8)
title('Amount of modification of the cross-section´s inertia')
ylabel('% of Gross cross-section´s Inertia')
xlabel('Load eccentricity (cm)')
legend('X axis')
hold on
plot(eccxyVec(:,2), IxyredVec(:,2)*100,'b -','LineWidth',1.8,...
    'DisplayName','Y axis')
hold on
% Adding labels
PuText=num2str(Pu);
PText=strcat('Pu = ',PuText,' Kg');
text(maxEc*0.8,max(IxyredVec(:,1))*100*0.7,PText)
AsText1=num2str(Aisr);
AsText=strcat('As = ',AsText1,' cm2');
text(maxEc*0.8,max(IxyredVec(:,1))*100*0.6,AsText)
%