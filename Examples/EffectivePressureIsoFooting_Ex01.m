% EffectivePressureIsoFooting_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the soil contact pressure's distribution of an isolated
%    RC footing. The function "RealPressuresFoot" is used.
%
%    To design the transversal footing's dimensions. The function
%    "designDimFootings" is used.
%
%    To design the footing's height by shear requirements. The function
%    "shearFootings" is used.
%
%    To determine the design bending moments at each of the footing's 
%    transversal cross-sections. The function "MomentDistributionFootings"
%    used.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-03-11
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clear all
clc

%% Basic material and geometry parameters
hefoot=25; % initial footing's height dimensions (cm)
rec=5; % concrete cover in all directions
d=hefoot-rec; % initial effective footing's height
bc=30; hc=30;  % transversal dimensions of the column
dimCol=[bc,hc]; 

fc=250; % concrete's compressive strneght f'c
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
     
%% Additional parameters
typeFooting=3; % Type of footing (1, 2 or 3) - see documentation
qadm=1.25; % Admissible bearing load of soil (Kg/cm2)
FS=1.2; % Safety Factor for the contact pressures
qu=qadm/FS; % Reduced soil's bearing load capacity
                         
%% Loads
load_cols=[1 -10.95e3 -11.8e5 -10.85e5]; % [n-load, Pu, Mux, Muy] Kg,cm
pu=load_cols(1,2);

%% Design of footing dimensions
% Note: in case it is required to design transversal dimensions
[be,le,contact_pressure]=designDimFootings(pu,qu,dimCol,hefoot,rec,...
                        typeFooting)

dim_zap=[be le];

%% Distribution of soil pressures
[qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_cols,be,le,typeFooting,...
                            dimCol,1);

%% Shear design
[d,qmax]=shearFootings(be,le,qprom,dimCol,pu,d,fc,typeFooting)
hefoot=d+rec; % Modified footing's height: minimum required by shear

%% Distribution of flexure over each footing's cross-section
dimpy=(be-hc-d)*0.5;
[mrL]=MomentDistributionFootings(qmax,dimpy,le)
dimpx=(le-bc-d)*0.5;
[mrB]=MomentDistributionFootings(qmax,dimpx,be)