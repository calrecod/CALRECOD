% PSO_Optim_RetainingRCWall_Ex01
%------------------------------------------------------------------------
% PURPOSE 
%    To optimally design a retaining reinforced concrete wall.
%
%------------------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

clc
clear all

%% Dimensions (geometry) - Optimization variables's range

H=220; % Total wall's height
D=0; % Back fill's depth

m1=1e15; % wall's stem's front slope
m2=1e15; % wall's stem's back slope

%% Materials

% Soil fill ----------------------------------------------------------
FiFill=30; % friction coefficient
wvFill=0.0018; % unit volume weight
beta=0; % surface slope

% Backing ground fill ------------------------------------------------
FiBackFill=30; % friction coefficient
alfa=0; % surface slope

% Foundation's soil
FiFound=30; % friction coefficient

% Concrete and reinforcement -----------------------------------------
fc=250; % concrete's compressive strength (Kg/cm2)
wvc=0.0024; % concrete unit volume weight
ductility=2; % ductility demand for the reinforcing steel
fy=4200; % yield stress of the reinforcing steel (Kg/cm2)

%% Design restrictions
qadm=1.25; % soil's bearing load capacity
minFSqadm=1.0; % safety factor against the soil's bearing load capacity

SlideSF=1.5; % safety factor against the sliding forces over the wall
TippingSF=2.0; % safety factor against tipping forces over the wall

% Rebar separation 
typeRebar=4; % eight-of-an-inch rebar diameter to be used 
sepMinRebars=10; % minimum rebar separation restriction to be complied
                 % for the reinforcement design of the wall

% Structural efficiency limit for each wall's element
maxEf=1.0;

%% Overload
qaf=0.01; % over the front fill (Kg/cm2)
qab=0; % over the back fill (Kg/cm2)

%% Support load
qs=34.67; % Kg/cm

%% Load Factors
LF_DL=1.3; % Dead Load Factor design

%% Optimization variables' search space
%      foot, heel, hf,  b
minDim=[20,  20,   15, 15];
MaxDim=[160, 160,  40, 40];

%% Optimization design function
[bestPerformance,bestPosition,besttippingFS,bestslideFS,...
bestLCap_FS,bestsepheel,bestefHeel,bestsepfoot,besteffoot,bestseptrunk,...
besteftrunk]=DesignRetainingWallPSO(minDim,MaxDim,H,D,m1,m2,FiFill,...
wvFill,beta,FiBackFill,alfa,FiFound,fc,fy,wvc,ductility,qadm,minFSqadm,...
SlideSF,TippingSF,typeRebar,sepMinRebars,maxEf,qaf,qab,qs,LF_DL)

%% Ploting optimal design
toe=bestPosition(1);
heel=bestPosition(2);
hf=bestPosition(3);
b=bestPosition(4);

plotRCWallDesign(H,m1,m2,toe,heel,hf,b,D,alfa,beta)