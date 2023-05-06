% Efficiency_AsymRebarDesign_Column_Ex_01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the structural efficiency of an assymmetrically
%    reinforced rectangular column cross-section in which only positive
%    bending moments are applied.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2023-03-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-04-16
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=60; %cm
h=60; %cm
concreteCover=[4 4]; % cm

%% Materials
fc=300; % Concrete compressive strength (Kg/cm2)
betac=0.85;

E=2.1e6; % Modulus of Elasticity of the reinforcing steel Kg/cm2
fdpc=fc*0.85; % (Kg/cm2)
fy=4200; % Yield stress of the reinforcing steel Kg/cm2

%% Additional parameters
npdiag=40; % number of points to be computed for the interaction diagrams

%% Loads
load_conditions=[1 -40000 38e5 53e5;
                 2 -21000 34e5 14e5]; % [nload, Pu, Mx, My] (Ton-m)

%% Rebar data
rebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];

numberRebars1=5;
numberRebars2=7;
numberRebars3=4;
numberRebars4=3;

% Total number of rebars placed over the cross-section
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4;

RebarTypeIndex1=5;
RebarTypeIndex2=4;
RebarTypeIndex3=6;
RebarTypeIndex4=7;

% Combination of rebar diameters (vector containing the 
% rebar diameters' indices for each of the four cross-section's boundary)
rebarcombo=[RebarTypeIndex1,RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4];
          
% Distribution of rebars over the cross-section
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarTypeIndex1,...
RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4);
      
%% Compute the column's interaction diagram
[diagrama,cPoints,Poc,Pot]=DiagramsAsymmetricRebar...
    (npdiag,rebarcombo,b,h,fy,fdpc,betac,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
    dispositionRebar);

%% Determine the column's structural resistance efficiency 
[maxef01,eficiencia01,cxy01]=effRecColsLinearSearch...
        (diagrama,load_conditions,Pot,Poc,cPoints) ;   
    
%% Plot the interaction diagram for better assessment 
rebarTypeslist=zeros(nv,1);
rebarTypeslist(1:numberRebars1)=RebarTypeIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarTypeIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarTypeIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarTypeIndex4;
  
diagramsFinalRebarCols(load_conditions,diagrama,dispositionRebar,...
                        h,b,rebarTypeslist);