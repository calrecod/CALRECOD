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

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clear all
clc

%% GEOMETRY
height=300; % column's length
b=30; % cross-section width
h=30; % cross-section height

%% MATERIALS
fy=4200; % yield stress of rebars 
fc=250; % concrete's compressive strength Kg/cm2
ws=0.0078; % volumetric weight of reinforcing steel (Kg/cm3)

%% Loads
load_conditions=[1 -53.52e3 -9.602e5 1e-5]; % [nload, Pu, Mx, My]
      
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

cols_sym_asym_isr="Symmetric";
if cols_sym_asym_isr=="Symmetric" || cols_sym_asym_isr=="Asymmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; % symmetric rebar

elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[28.93];
end

condition_cracking="Cracked"; % "Cracked" or "Non-cracked"

% Plot options: 
optimaPlot=1; % optima reinforcement area convergence plot
plotsISRdiagrams=1; % interaction diagrams
plotRebarDesign=1; % detailed reinforced cross-section

% Ductility demand
ductility=3;
             
%% Optimization process
[Inertia_xy_modif,b,h,bestArrangement,bestdisposicionRebar,cost_elem_col,...
Ac_sec_elem,Ef_sec_col,Mr_col]=isrColumnsSymAsym(pu_cols,height,b,h,rec,...
fy,fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,...
ws,RebarAvailable,optimaPlot,plotsISRdiagrams,plotRebarDesign);

%% Export results

dimColumnsCollection=[b h height rec(1)];
nbarColumnsCollection=length(bestdisposicionRebar(:,1));

directionData='C:\Users\lfver\OneDrive\Desktop\OneDrive\DynamoRevit\Dynamo_visualization\';

ExportResultsColumn(directionData,dimColumnsCollection,bestdisposicionRebar,...
    nbarColumnsCollection,bestArrangement,cols_sym_asym_isr,[0,0,0],1);
