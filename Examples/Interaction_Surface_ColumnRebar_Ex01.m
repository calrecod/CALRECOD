% Interaction_Surface_ColumnRebar_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the interaction surface of a rebar reinforced 
%    rectangular column.
%    
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=60; % cross-section width (cm)
h=60; % cross-section height (cm)

% Materials
fc=300; % concrete's compressive strength kg/cm2

E=2.1e6; % Modulus of Elasticity of the reinforcing steel (kg/cm2)

fdpc=fc*0.85; % reduces f'c (Kg/cm2)
fy=4200; % Yield stress of the reinforcing steel (Kg/cm2)
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

%% Additional parameters
npdiag=40;
concreteCover=[4 4]; %cm

%% Rebar data
% Database of the commercially available rebar
rebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];

% Number of rebars at of the four cross-section's boundaries
numberRebars1=5;
numberRebars2=7;
numberRebars3=4;
numberRebars4=3;

% Calculation of rebar area at each cross-section boundary
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4;

% Rebar diameter indices at each of the four cross-section's boundaries
RebarIndex1=5;
RebarIndex2=4;
RebarIndex3=6;
RebarIndex4=7;

% eight-of-an-inch (or rebar diameter) rebar at each cross-section boundary
comborebar=[RebarIndex1,RebarIndex2,RebarIndex3,RebarIndex4];

% Compute the distribution of rebars over the cross-section (local rebar 
% coordinates) 
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarIndex1,...
RebarIndex2,RebarIndex3,RebarIndex4);
     
%% Interaction surface
[supX,supY,supZ]=InteracSurfaceColRec(b,h,comborebar,npdiag,fy,fdpc,...
beta1,E,numberRebars1,numberRebars2,numberRebars3,numberRebars4,...
rebarAvailable,dispositionRebar);
% 
