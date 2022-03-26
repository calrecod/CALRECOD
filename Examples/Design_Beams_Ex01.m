clc
clear all

pu_beams=38.6; % unit construction assembly cost of steel reinforcement
duct=3; % high ductility demand
span=500; % cm
b=30; % width (cm)
h=60; % height (cm)
                
h_rec_sections=[5 5 3 3 5 5]; % [rec_left_up, rec_left_low, rec_mid_up,
                              % rec_mid_low, rec_right_up, rec_right_low]
                              
b_rec=3; % lateral concrete cover
fc=280; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
load_conditions=[1 -33.0 29.0 -31.0]; %Ton-m
cols_sym_asym_isr="Symmetric";
plots=1;

[sep_bars,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
    barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
    barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
    barArrangementIzqComp,minAreaVar_3sec,Ef_elem_sec_t,bestCostVar,ef_var,...
    minAreaVar_prom,Mr_3section]=beamsISR(pu_beams,span,b,h,h_rec_sections,...
    fc,fy,load_conditions,cols_sym_asym_isr,duct,b_rec,plots)
