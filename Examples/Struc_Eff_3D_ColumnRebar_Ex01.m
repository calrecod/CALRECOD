% Struc_Eff_3D_ColumnRebar_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the structural efficiency of a biaxially loaded rectangular
%    column cross-section through a 3D analysis with its interaction
%    surface. 
%    
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------


clc
clear all

%% Geometry
b=60; %cm
h=60; %cm

%% Materials
fc=300; %kg/cm2
beta1=0.85;
fdpc=fc*0.85; % (Kg/cm2)
fy=4200; % Kg/cm2
if fc<280
    beta1=0.85;
elseif fc>=280
    beta1=1.05-fc/1400;
    if (beta1<0.65)
        beta1=0.65;
    elseif (beta1>0.85)
        beta1=0.85;
    end
end
E=2.1e6; % Modulus of Elasticity of the reinforcing steel (kg/cm2)

%% Additional parameters
npdiag=40; % Number of points of the interaction diagram
concreteCover=[4 4]; % cm

%% Loads
load_conditions=[1 -15000 -55e5 32e5;
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
numberRebars2=7;
numberRebars3=4;
numberRebars4=3;

% Total numbers of rebars placed over the cross-section
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4;

RebarTypeIndex1=5;
RebarTypeIndex2=4;
RebarTypeIndex3=6;
RebarTypeIndex4=7;

% Rebar diameter over each of the four cross-section's boundaries
comborebar=[RebarTypeIndex1,RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4];

% Compute the distribution of rebars over the cross-section (local rebar 
% coordinates) 
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarTypeIndex1,...
RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4);
            
%% Interaction diagram in the load directions (3D Interaction surfaces)
[tablaEff01,iloadmax,trans_load_condition,gamma,diagramIntAxis1,...
    newdispositionRebar,section,cmax,CP]=multiDiagAxisColRec(b,h,...
    load_conditions,comborebar,npdiag,fy,fdpc,beta1,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
    dispositionRebar);

% Plotting results
rebarTypeslist(1:numberRebars1)=RebarTypeIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarTypeIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarTypeIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarTypeIndex4;
  
PlotRotRecSecRebarCols(section,trans_load_condition,diagramIntAxis1,...
                    newdispositionRebar,b,h,rebarTypeslist);
