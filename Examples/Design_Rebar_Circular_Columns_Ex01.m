% Design_Rebar_Circular_Columns_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to savings in reinforcing area)
%    a circular column' cross-section (considering uniaxial bending)
%    only with respect to one direction of action (left or right) for 
%    the load combinations - two options of reinforcement 
%    optimization are possible: with a pure ISR, symmtrical rebar
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-01-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% GEOMETRY  
height=400; %column's length (cm)
diam=50; % cross-section diameter (cm)
rec=4; % concrete cover (cm)

%% MATERIAL
fy=4200; % yield stress of rebars 
E=fy/0.0021; % Modulus of Elasticity of the reinforcing steel (Kg/cm2)
fc=280; % concrete's compressive strength (Kg/cm2)
wac=7.8e-3; % unit volume weight of reinforcing steel (Kg/cm2)

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
load_conditions=[1 -15e3 32e5]; % [nload, Pu, Mu] (Kg-cm)

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
if cols_rebar_isr=="Rebar"
    
    plotRebarDesign=1; % To visualize the rebar design layout plot
    npdiag=40;
    [MrColRebar,Inertia,bestArea,bestCost,bestdiagram,bestnv,bestobar,...
    bestEf,bestArrangement,bestDisposition,bestc]=optimalRebarCirc...
    (diam,rec,act,E,npdiag,0.85*fc,puColsRebar,load_conditions,...
    barsAvailable,wac,height,plotRebarDesign);
end

%% EXPORT RESULTS
coordBaseCols=[0 0 0]; % global location of the column base [x,y,z]
dimColumnsCollection=[diam height rec];
nbarColumnsCollection=length(bestDisposition(:,1));
directionData='C:\Users\luizv\OneDrive\DynamoRevit\Dynamo_visualization\';

ExportResultsColumnCirc(directionData,dimColumnsCollection,bestDisposition,...
    nbarColumnsCollection,bestArrangement,cols_rebar_isr,coordBaseCols)
