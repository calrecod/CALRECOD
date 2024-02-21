% StabilitySlenderRCols_Rebar_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To analyse the variation in structural efficiency of an asymmetrically
%    reinforced rectangular cross-section when considering 
%    cracking-cross-section mechanisms and Geometrical Non-linearities
%
%----------------------------------------------------------------

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
k=1.5; % slederness factor
npdiag=50; % number of points for the computation of the interaction
           % diagram

concreteCover=[4 4]; %cm

%% Loads
load_conditions=[1 -15000 28e5 22e5]; % [nload, Pu, Mx, My]
P=load_conditions(2);
Mx=load_conditions(3);
My=load_conditions(4);
Vx=0.0; % Lateral loads
Vy=0.0;

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
numberRebars1=2;
numberRebars2=8;
numberRebars3=0;
numberRebars4=0;

% Total number of rebars placed at each cross-section boundary
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4;

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
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarIndex1,...
RebarIndex2,RebarIndex3,RebarIndex4);
   
%% Computation of the Plastic Center with respect to the X-axis
rebarTypeslist(1:numberRebars1)=RebarIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarIndex4;
 
[PCX]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)
                        
%% Computation of the Plastic Center with respect to the Y-axis
rebarCoordinates=[dispositionRebar(:,1) dispositionRebar(:,2)]; 

% Invert rebar local coordinates (the cross-section is rotated 90°)
dispositionRebar(:,1)=-rebarCoordinates(:,2);
dispositionRebar(:,2)=rebarCoordinates(:,1);

% Invert cross-section dimensions
h=dimensionsColumn(1);
b=dimensionsColumn(2);  

% Invert rebar diameters over the cross-section
rebarTypeslist2(1:numberRebars3)=RebarIndex3;
rebarTypeslist2(numberRebars3+1:numberRebars3+numberRebars4)=RebarIndex4;
rebarTypeslist2(numberRebars3+numberRebars4+1:numberRebars3+numberRebars4+...
    numberRebars2)=RebarIndex2;
rebarTypeslist2(numberRebars3+numberRebars4+numberRebars2+1:nv)=RebarIndex1;
 
[PCY]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist2,...
                            rebarAvailable)   

%% Interaction diagrams with respect to the X and Y axis
h=dimensionsColumn(2);
b=dimensionsColumn(1); 

[diagrama,cPoints,Poc,Pot]=DiagramsAsymmetricRebar...
    (npdiag,rebarcombo,b,h,fy,fdpc,betac,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
    rebarCoordinates);

%% Structural resistance efficiency according to applied load combinations
[maxef01,eficiencia01,cxy01]=effRecColsLinearSearch...
        (diagrama,load_conditions,Pot,Poc,cPoints)     

%% Plot the interaction diagram for better assessment 
diagramsFinalRebarCols(load_conditions,diagrama,rebarCoordinates,...
                        h,b,rebarTypeslist);

%% Modified inertia momentum of the cross-section

% Cracking mechanisms are considered with the "transformed 
% cross-section method"
eccentricityXY=[abs(Mx/P),abs(My/P)]; % loads eccentricities
conditionCrack="Cracked";
PC=[PCX,PCY];

% Total rebar area at each cross-section's boundary
ab1=rebarAvailable(RebarIndex1,2)^2*pi/4*numberRebars1;
ab2=rebarAvailable(RebarIndex2,2)^2*pi/4*numberRebars2;
ab3=rebarAvailable(RebarIndex3,2)^2*pi/4*numberRebars3;
ab4=rebarAvailable(RebarIndex4,2)^2*pi/4*numberRebars4;

t1bar=ab1/(b-2*concreteCover(1));
t2bar=ab2/(b-2*concreteCover(1));
t3bar=ab3/(h-2*concreteCover(2));
t4bar=ab4/(h-2*concreteCover(2));

InertiaXY=[1/12*b*h^3,1/12*h*b^3]; % gross cross-section's inertia

[InertiaXY,Atransfxy,elimxy]=CrackingColumnsAsym(h,b,fdpc,concreteCover,...
    eccentricityXY,t1bar,t2bar,t3bar,t4bar,P,cxy01,conditionCrack,PC);

%% Euler's critical load: 
% the smaller inertia value of both axis directions must be taken as reference
if InertiaXY(1)<=InertiaXY(2)
    I=InertiaXY(1);
else
    I=InertiaXY(2);
end

% Critical Euler's load
Pcr=pi^2*Ec*I/(k*height)^2;

%% Amplification bending moments for each axis direction
[Deltax,Mampx]=MomAmpColsGeomNL(fc,k,InertiaXY(1),height,Vx,P,...
                                Mx,b,h,0); % In the x-axis

h=dimensionsColumn(1);
b=dimensionsColumn(2);
[Deltay,Mampy]=MomAmpColsGeomNL(fc,k,InertiaXY(2),height,Vy,P,...
                                My,b,h,0); % In the y-axis

%% Amplified load bending moments
load_conditions_modif=[1 P Mampx Mampy]

%% Reduced resistance efficiency due to the amplified moments
% by considering geometrical non-linearities (P-Delta effects)
[maxef02,eficiencia02,cxy02]=effRecColsLinearSearch...
        (diagrama,load_conditions_modif,Pot,Poc,cPoints)

diagramsFinalRebarCols(load_conditions_modif,diagrama,rebarCoordinates,...
                        h,b,rebarTypeslist);