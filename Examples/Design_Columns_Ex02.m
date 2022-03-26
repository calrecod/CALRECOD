clear all
clc

cols_sym_asym_isr="Asymmetric"; % Distribution of reinforcement
if cols_sym_asym_isr=="Symmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; % symmetric rebar
elseif cols_sym_asym_isr=="Asymmetric"
    pu_cols=[32.23, 32.10, 31.96, 31.96, 31.96, 31.96, 31.96;
            36.0, 36.0, 36.0, 36.0, 36.0, 36.0, 36.0]; % asymmetric rebar

elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[28.93];
end

height=400; % column length or height (cm)
b=60; % width cross-section (cm)
h=60; % height cross-section (cm)
fy=4200; % yield stress of rebars 
fc=280; % Kg/cm2
load_conditions=[1 15 32 8]; % [nload, Pu, Mx, My] (Ton-m)
rec=[4 4]; % concrete cover: [coverx covery] (cm)

condition_cracking="Cracked";
ductility=3; % high ducitlity demand on cross-section

optimaPlot=1;
plotISRdiagrams=1;
plotRebarDesign=1;
% _______________________________________________________________________

[Inertia_xy_modif,b,h,bestArrangement,best_disposicion,cost_elem_col,...
Ac_sec_elem,Ef_sec_col,Mr_col]=isrColumnsSymAsym(pu_cols,height,b,h,rec,...
fy,fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,optimaPlot,...
plotISRdiagrams,plotRebarDesign);