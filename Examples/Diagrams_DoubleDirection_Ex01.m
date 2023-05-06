% Diagrams_DoubleDirection_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the interaction diagrams of an asymmetrically rein-
%    forced rectangular column's section for both load-directions
%    with respect of each cross-section local axis
%
%    Note: function dispositionRebarAsymmetric is the one required
%          to compute the asymmetrical distribution of rebars
%          (local rebar coordinates over the cross-section)
%
%          function EvalAsymDoubleDirection is the noe required to
%                   compute the interaction diagrams
%
%          function effRecColsDoubleDirecBS is used to assess the
%                   resistance efficiency of the cross-section against
%                   the load conditions
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-08-02
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

%% Geometry
b=60; %cm
h=60; %cm

%% Materials
fc=300; % kg/cm2
betac=0.85;

E=2.1e6; % kg/cm2

fdpc=fc*0.85; % (Kg/cm2)
fy=4200; % Kg/cm2

%% Additional parameters
concreteCover=[4 4]; % cm
npdiag=40; % number of points to compute for the interaction diagrams

%% Loads
load_conditions=[1 -15000 -25e5 22e5;
                 2 -40000 38e5 53e5;
                 3 -31000 12e5 -5.5e5]; % [nload, Pu, Mx, My] 
             
%% Rebar data
% Database of the commercially available rebar
rebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];

numberRebars1=5;
numberRebars2=2;
numberRebars3=6;
numberRebars4=1;

% Total number of rebars placed over the cross-section
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4;

RebarIndex1=5;
RebarIndex2=4;
RebarIndex3=6;
RebarIndex4=7;

rebarcombo=[RebarIndex1,RebarIndex2,RebarIndex3,RebarIndex4];
         
% Distribution of rebars over the cross-section (local rebar 
% coordinates) 
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarIndex1,...
RebarIndex2,RebarIndex3,RebarIndex4);

%% Compute the column's interaction diagram
[diagramaInteraccion1,diagramaInteraccion2,pot,poc,cp_axis,...
cvector1,cvector2]=EvalAsymDoubleDirection(npdiag,rebarcombo,b,h,...
fy,fdpc,betac,E,numberRebars1,numberRebars2,numberRebars3,...
numberRebars4,rebarAvailable,dispositionRebar,concreteCover);

%% Determine the column's structural efficiency
[maxef,tablaEficiencias,cxy]=effRecColsDoubleDirecLS(diagramaInteraccion1,...
           diagramaInteraccion2,load_conditions,pot,poc,cvector1,cvector2)
       
%% Plot the interaction diagram for better assessment 

% A list of the rebar diameters' index is created. Size=n-rebars x 1, 
% starting with the rebars at the top of the cross-section, then at the 
% bottom, then at the left side and finally at the right side
rebarTypeslist(1:numberRebars1)=RebarIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarIndex4;

diagDoubleDirecAsymRebarCols(load_conditions,diagramaInteraccion1,...
    diagramaInteraccion2,dispositionRebar,h,b,rebarTypeslist);
                    