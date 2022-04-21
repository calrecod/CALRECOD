clear all
clc

fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
qadm=2.5; % Admissible bearing load of soil (Kg/cm2)
FS=1.5; % Design Safety Factor
qu=qadm*FS; % Design bearing load of soil
pu_steel_footings=26.75; % unit construction cost of reinforcement assembly
     
he_footing=40; % initial footing height dimensions (cm)
fc=300; %f'c
rec=5; % concrete cover in all directions
bc=30; hc=50;  % transversal dimensions
dimCol=[bc,hc]; % supporting column's dimensions

load_conditions=[1 20 35 24]; % Ton,m
pu=load_conditions(1,2);

% Design of footing dimensions___________________________________________
% Note: in case it is required to design transversal dimensions
[be,le,contact_pressure]=designDimFootings(pu,qu,dimCol,he_footing,rec);

dim_zap=[be le];

cols_sym_asym_isr="Rebar"; % ''ISR'' or anything else (''Symmetric'', ''Rebar'')
ductility=3; % to set the MAX reinforcement quantity (1, 2 or 3)

optimConv=1;
PlotRebarDesign=1;

[hmodif,m_max_eje,barDispositionFootings,arrangement_bar_footings,...
nbars_footings,AcBar,bestCost_elem,list_ef_footings,list_mr_footings]=...
isrFootings(pu_steel_footings,he_footing,dim_zap(1),dim_zap(2),...
rec,fc,fy,load_conditions,dimCol,cols_sym_asym_isr,ductility,optimConv,PlotRebarDesign)
