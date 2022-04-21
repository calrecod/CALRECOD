clear all
clc

% GEOMETRY __________________________________________________________
height=400; %cm
b=60; % cm
h=60; % cm
fy=4200; % yield stress of rebars 
fc=280; % Kg/cm2
load_conditions=[1 14 22 31]; % [nload, Pu, Mx, My] (Ton-m)
rec=[4 4]; % [coverx covery] (cm)

% ADITIONAL PARAMETERS ______________________________________________
cols_sym_asym_isr="Symmetric";
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
