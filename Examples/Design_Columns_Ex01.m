% Design_Columns_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing area)
%    a rectangular column' cross-section, considering biaxial bending
%    loads - three options of reinforcement optimization are
%    possible: with a pure ISR, symmtrical rebar or asymmetrical rebar
%
%    Note: function "isrColumnsSymAsym" is the only on required.
%
%          function "ExportResultsColumn" is the one required to export the
%          results to given disc direction through the parameter 
%          "directionData" as a string variable
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clear all
clc

%% GEOMETRY
height=300; % column's length
b=30; % cross-section width
h=40; % cross-section height

%% MATERIALS
fy=4200; % yield stress of rebars 
fc=250; % concrete's compressive strength Kg/cm2
ws=0.0078; % volumetric weight of reinforcing steel (Kg/cm3)

%% Loads
load_conditions=[1 -53.52e3 -9.602e5 10e5]; % [nload, Pu, Mx, My]
      
%% Rebar data
% Commerically available rebars
RebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
%% ADITIONAL PARAMETERS
rec=[4 4]; % [coverx covery] - concrete cover

cols_sym_asym_isr="Asymmetric";

pu_cols=[28.93]./18;
puColBuild=[1.2,0.9,128,214.75,0.04,0.105,7];

condition_cracking="Cracked"; % "Cracked" or "Non-cracked"
dataCFA=[0,1,1,2]; % Data required for the computation of the 
                   % constructability factor for the rebar designs of columns
                   % Lower and upper values for the range of the CFA and 
                   % weight factors for the uniformity of rebars and rebar 
                   % diameters
% Plot options: 
optimaPlot=1; % optima reinforcement area convergence plot
plotsISRdiagrams=1; % interaction diagrams
plotRebarDesign=1; % detailed reinforced cross-section

% Ductility demand
ductility=3;
             
%% Optimization process
[Inertia_xy_modif,b,h,bestArrangement,bestdisposicionRebar,cost_elem_col,...
Ac_sec_elem,Ef_sec_col,Mr_col,bestCFA]=isrColumnsSymAsym(pu_cols,height,b,h,rec,...
fy,fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,...
ws,RebarAvailable,puColBuild,dataCFA,optimaPlot,plotsISRdiagrams,plotRebarDesign);

%% Export results

dimColumnsCollection=[b h height rec(1)];
nbarColumnsCollection=length(bestdisposicionRebar(:,1));

% disc's route to export design results
directionData='C:\Users\lfver\OneDrive\Desktop\OneDrive\DynamoRevit\Dynamo_visualization\';

ExportResultsColumn(directionData,dimColumnsCollection,bestdisposicionRebar,...
    nbarColumnsCollection,bestArrangement,cols_sym_asym_isr,[0,0,0]);
