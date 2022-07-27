% Design_Columns_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing area)
%    a rectangular column' cross-section (considering biaxial bending)
%    only with respect to one direction of action for the load combinations
%    (left or right) - three options of reinforcement optimization are
%    possible: with a pure ISR, symmtrical rebar or asymmetrical rebar
%
%    Note: isrColumnsSymAsym function is the only on required.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clear all
clc

% GEOMETRY __________________________________________________________
height=400; %cm
b=60; % cm
h=60; % cm
fy=4200; % yield stress of rebars 
fc=280; % Kg/cm2
load_conditions=[1 15 32 8]; % [nload, Pu, Mx, My] (Ton-m)
rec=[4 4]; % [coverx covery] (cm)

% ADITIONAL PARAMETERS ______________________________________________
cols_sym_asym_isr="ISR";
if cols_sym_asym_isr=="Symmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; % symmetric rebar
elseif cols_sym_asym_isr=="Asymmetric"
    pu_cols=[32.23, 32.10, 31.96, 31.96, 31.96, 31.96, 31.96;
            36.0, 36.0, 36.0, 36.0, 36.0, 36.0, 36.0]; % asymmetric rebar

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

% _______________________________________________________________________

[Inertia_xy_modif,b,h,bestArrangement,best_disposicion,cost_elem_col,...
Ac_sec_elem,Ef_sec_col,Mr_col]=isrColumnsSymAsym(pu_cols,height,b,h,rec,...
fy,fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,...
optimaPlot,plotsISRdiagrams,plotRebarDesign);
