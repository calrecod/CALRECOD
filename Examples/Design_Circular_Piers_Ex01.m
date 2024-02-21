% Design_Circular_Piers_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to savings in reinforcing area)
%    a circular column' cross-section. Rebar designs in packages of two are
%    considered.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-01-24
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% GEOMETRY  
height=500; %column's length (m)
diam=150; % cross-section diameter (m)
rec=4; % concrete cover (cm)

%% MATERIAL
fy=4200; % yield stress of rebars 
Es=fy/0.0021; % Modulus of Elasticity of the reinforcing steel (Kg/cm2)
fc=500; % concrete's compressive strength (Kg/cm2)
wac=7.8e-3; % unit volume weight of reinforcing steel (Kg/cm2)
esy=fy/(1.15*Es);
Ec=12000*sqrt(fc);

%% ADITIONAL PARAMETERS 
cols_rebar_isr="Rebar";
puColsISR=[28.93];

% symmetrical rebar
puColsRebar=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; 

% Plot options: 
optimaConvPlot=1; % optima reinforcement area convergence plot
plotsISRdiagrams=1; % interaction diagrams

% Ductility demand
ductility=3;

%% Loads
load_conditions=[1 -7.159e3/9.81 sqrt((4576/9.81*100*1000)^2+(1270/9.81*100*1000)^2);
                 2 -7.500e3/9.81 sqrt((3720/9.81*100*1000)^2+(1296/9.81*100*1000)^2);
                 3 -7.238e3/9.81 sqrt((713/9.81*100*1000)^2+(4355/9.81*100*1000)^2);
                 4 -7.082e3/9.81 sqrt((456/9.81*100*1000)^2+(4355/9.81*100*1000)^2)]; % [nload, Pu, Mu] (Kg-cm)

%% Rebar data
% Commerically available rebars
barsAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
%% ISR OPTIMIZATION
[cost_elem_col,act,Ef_sec_col,MrtCol,t_value,c]=isrCircCols(puColsISR,...
height,wac,diam,rec,fy,fc,load_conditions,ductility,optimaConvPlot,...
plotsISRdiagrams);

%% REBAR OPTIMIZATION (OPTIONAL)

plotRebarDesign=1; % To visualize the rebar design layout plot
npdiag=40;
[MrColRebar,Inertia,bestArea,bestCost,bestdiagram,bestnv,bestobar,...
bestEf,bestArrangement,bestDisposition,bestc,criticalEf]=optimalRebarCirc2pack...
(diam,rec,act,Es,npdiag,0.85*fc,puColsRebar,load_conditions,...
barsAvailable,wac,height,plotRebarDesign);

                                         
