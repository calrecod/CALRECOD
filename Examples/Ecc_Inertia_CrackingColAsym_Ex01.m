% Ecc_Inertia_CrackingColAsym_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the modification of the inertia of a rectangular column
%    cross-section reinforced with asymmetrical rebar subject to biaxial 
%    bending-compression forces as the load eccentricities increase
%
%    Note: effRecColsLinearSearch function is used to obtain the interaction
%          diagrams with respect to both cross-section axis directions, as
%          well as the neutral axis depth for each direction according to
%          the load combinations
%          
%          CrackingColumnsAsym function is the one used to compute the
%          modified inertia momentums of the cross-section for axis
%          directions, according the laod combinations (load eccentricity)
%          and neutral axis depth c
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=40;
h=40;
dimensionsColumn=[b h];
height=500; % Column's length

%% Materials
fc=300;
Ec=14000*sqrt(fc); % Modulus of Elasticity of the concrete

betac=0.85;
E=2.1e6;%kg/cm2

fdpc=fc*0.85; % reduced concrete's compressive strength (Kg/cm2)
fy=4200; % Yield stress of reinforcing steel Kg/cm2

%% Additional structural parameters
k=2; % slederness factor
npdiag=90; % number of points for the computation of the interaction
           % diagram

concreteCover=[4 4]; %cm

%% Rebar data

% Commerically available rebar diameters (eigh-of-an-inch)
rebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];
                
% Number of rebars at each cross-section boundary
nb1=2;
nb2=2;
nb3=0;
nb4=6;

% Total number of rebars placed at each cross-section boundary
nv=nb1+nb2+nb3+nb4;

% Rebar diameter's indices (from the "rebarAvailable" array) at each 
% cross-section boundary
RebarIndex1=4;
RebarIndex2=4;
RebarIndex3=4;
RebarIndex4=4;
  
rebarcombo=[RebarIndex1,RebarIndex2,RebarIndex3,RebarIndex4];

% Rebar coordinates
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,nb1,nb2,nb3,nb4,rebarAvailable,RebarIndex1,...
RebarIndex2,RebarIndex3,RebarIndex4);
   
%% Computation of the Plastic Center with respect to the X-axis
rebarTypeslist(1:nb1)=RebarIndex1;
rebarTypeslist(nb1+1:nb1+nb2)=RebarIndex2;
rebarTypeslist(nb1+nb2+1:nb1+nb2+nb3)=RebarIndex3;
rebarTypeslist(nb1+nb2+nb3+1:nv)=RebarIndex4;
 
[PCX]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable);
                        
%% Computation of the Plastic Center with respect to the Y-axis
rebarCoordinates=[dispositionRebar(:,1) dispositionRebar(:,2)]; 

% Invert rebar local coordinates (the cross-section is rotated 90°)
dispositionRebar(:,1)=-rebarCoordinates(:,2);
dispositionRebar(:,2)=rebarCoordinates(:,1);

% Invert cross-section dimensions
h=dimensionsColumn(1);
b=dimensionsColumn(2);  

% Invert rebar diameters over the cross-section
rebarTypeslist2(1:nb3)=RebarIndex3;
rebarTypeslist2(nb3+1:nb3+nb4)=RebarIndex4;
rebarTypeslist2(nb3+nb4+1:nb3+nb4+...
    nb2)=RebarIndex2;
rebarTypeslist2(nb3+nb4+nb2+1:nv)=RebarIndex1;
 
[PCY]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist2,...
                            rebarAvailable);
                        
PC=[PCX,PCY];
h=dimensionsColumn(2);
b=dimensionsColumn(1);

%% Interaction diagrams with respect to the X and Y axis
[diagrama,cPoints,Poc,Pot]=DiagramsAsymmetricRebar...
    (npdiag,rebarcombo,b,h,fy,fdpc,betac,E,nb1,...
    nb2,nb3,nb4,rebarAvailable,rebarCoordinates);

%% Additional computations
a1bar=nb1*rebarAvailable(RebarIndex1,2)^2*pi/4;
a2bar=nb2*rebarAvailable(RebarIndex2,2)^2*pi/4;
a3bar=nb3*rebarAvailable(RebarIndex3,2)^2*pi/4;
a4bar=nb4*rebarAvailable(RebarIndex4,2)^2*pi/4;

t1bar=a1bar/(b-2*concreteCover(1));
t2bar=a2bar/(b-2*concreteCover(1));
t3bar=a3bar/(h-2*concreteCover(2));
t4bar=a4bar/(h-2*concreteCover(2));

Abar=a1bar+a2bar+a3bar+a4bar % proposed reinforcement area
Amin=0.01*b*h % Min reinforcement area by code
Amax=0.04*b*h % Max reinforcement area by code

Inertia0=[b*h^3/12,h*b^3/12]; % gross section's inertia

%% Loads
Pu=-45.3e3; % Compression load (Kg)

%% Additional parameters
conditionCrack="Cracked";

%% Main loop
IxyredVec=[];
eccxyVec=[];
maxEc=60; % Max load eccentricity for both section's axis
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
    [maxef01,eficiencia01,cxy01]=effRecColsLinearSearch...
        (diagrama,conditions,Pot,Poc,cPoints);
   
    %% Modified momentum of inertia for both axis directions
    % cracking mechanisms are considered through the 
    % "transform section" method

    [InertiaXYmodif,Atransfxy,elimxy]=CrackingColumnsAsym(h,b,fc*0.85,...
    concreteCover,eccentricityXY,t1bar,t2bar,t3bar,t4bar,Pu,cxy01,...
    conditionCrack,PC);

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
text(maxEc*0.8,max(IxyredVec(:,1))*100*0.9,PText)
AsText1=num2str(Abar);
AsText=strcat('As = ',AsText1,' cm2');
text(maxEc*0.8,max(IxyredVec(:,1))*100*0.7,AsText)
%