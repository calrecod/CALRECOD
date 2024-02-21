% Design_Retaining_RC_Wall_Ex01
%------------------------------------------------------------------------
% PURPOSE 
%    To analyse and design the reinforcement of reinforced concrete
%    retaining wall for a given set of fixed dimensions and soil's
%    mechanical properties.
%
%------------------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

clc
clear all

%% Dimensions (geometry)

H=220; % total wall's height
D=0; % depth of the back soil

m1=1e15; % wall's front stem's slope
m2=1e15; % wall's back stem's slope

toe=25; % toe's length
heel=120; % heel's length
hf=25; % width of toe and heel
b=20; % stem's upper width

%% Materials

% Concrete and reinforcement -----------------------------------------
fc=250;
wvc=0.0024; % unit volume weight of the concrete
ductility=2; % ductility demand of the reinforcing steel

fy=4200; % yield stress of the reinforcing steel

% Soil fill ----------------------------------------------------------
FiFill=30; % friction coefficient of the front soil fill
wvFill=0.0018; % unit volume weight of the front soil fill
beta=10; % fill's top slope

% Backing ground fill ------------------------------------------------
FiBackFill=30; % friction coefficient of the back soil fill
alfa=0; % back fill's top slope

% Foundation's soil --------------------------------------------------
FiFound=30; % friction coefficient of the foundation soil

%% Design restrictions
qadm=1.25; % soil's bearing load capacity (Kg/cm2)
minFSqadm=1.0; % limit of the soil's load capacity safety factor  

SlideSF=1.1; % safety factor against sliding forces over the wall
TippingSF=2.0; % safety factor against tipping forces over the wall

% Rebar separation ---------------------------------------------------
typeRebar=4; 
sepMinRebars=10; % minimum rebar separation restriction for the 
                 % reinforcement on each wall's element

% Structural efficiency ----------------------------------------------
maxEf=1.0; % limit of the structural efficiency of each 
           % reinforced wall's element

%% Overload
qaf=0.01; % over the front soil's fill
qab=0; % over the back soil's fill

%% Support load
qs=34.67; % vertical linear load acting along the top of the wall's stem

%% Load Factors
LF_DL=1.3; % Dead Load Factor for design

%% Design function
[compliedRestric,areaWall,linearWeigthWall,tippingFS,slideFS,LCap_FS,...
sepheel,efHeel,sepfoot,effoot,septrunk,eftrunk]=RetainingRCWall(toe,...
heel,hf,b,FiFill,H,D,m1,m2,wvFill,beta,FiBackFill,alfa,FiFound,fc,fy,wvc,...
ductility,qadm,minFSqadm,SlideSF,TippingSF,typeRebar,sepMinRebars,maxEf,...
qaf,qab,qs,LF_DL)

%% Ploting results
plotRCWallDesign(H,m1,m2,toe,heel,hf,b,D,alfa,beta)