% Optimal_Design_RecBeams_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing area)
%    a beam element for all its three critical cross-sctions (left,middle
%    right)
%
%    Note: beamsISR function is the only on required.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-04-30
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

%% Geometry 
span=500; % cm
b=30; % width (cm)
h=60; % height (cm)
                
h_rec=5; % 
b_rec=3; % lateral concrete cover

%% Materials
fc=280; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
wac=7.8e-3; % unit volume weight of the reinforcing steel (Kg/cm3)

%% Load conditions
load_conditions=[1 -33.0e5 29.0e5 -31.0e5]; %Kg-cm (flexure)
shear_beam=linspace(12e3,-22e3,7); %Kg (shear)

%% Additional parameters
cols_sym_asym_isr="Rebar"; % ''Rebar'' or ''ISR''

% Plot options
plots=1; % for reinforced cross-section plotts (0-No,1-Yes)
graphConvergence=1; % for optimal ISR area convergence (0-No,1-Yes)

pu_beams=38.6; % unit construction assembly cost of steel reinforcement
duct=1; % high ductility demand

%% Rebar data
% Available commercial rebar diameters (in eight-of-an-inch)
                %type diam
rebarAvailable=[4 4/8.*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];

%% OPTIMAL DESIGN 
[sep_bars,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
barArrangementIzqComp,minAreaVar_3sec,Ef_elem_sec_t,bestCostVar,ef_var,...
minAreaVar_prom,Mr_3section]=beamsISR(pu_beams,span,wac,b,h,h_rec,...
rebarAvailable,fc,fy,load_conditions,cols_sym_asym_isr,duct,b_rec,plots,...
graphConvergence);

%% SHEAR DESIGN
rho=(sum(minAreaVar_3sec)/length(minAreaVar_3sec))/(b*h);
[s1,s2,s3,d1,d3]=shearDesignBeams(span,b,h,h_rec,rho,fc,fy,shear_beam);
                     
%% Exporting results
dim_beams_collection(1,:)=[b h span h_rec]; % general geometry data
shear_beam_design_collec=[s1,s2,s3,d1,d3]; % shear design data

disposition_rebar_beams3sec=[dispositionBar_Izq; % rebar design coordinates
                             dispositionBar_Center;
                             dispositionBar_Der];
                         
              % number of rebar for each cross-section
              % both in tension and compression
nbar_beams_collection(1,:)=[length(barArrangementIzqTens)...
                            length(barArrangementIzqComp)...
                            length(barArrangementCentralTens)...
                            length(barArrangementCentralComp)...
                            length(barArrangementDerTens)...
                            length(barArrangementDerComp)]; 
                        
              % rebar diameters at each cross-section
              % in tension and compression
beamsDiamIndex=[barArrangementIzqTens;
                    barArrangementIzqComp;
                    barArrangementCentralTens;
                    barArrangementCentralComp;
                    barArrangementDerTens;
                    barArrangementDerComp];
                    
collection_beams_diamIndex=[beamsDiamIndex];
                        
coordEndBeams=[0 0 300]; % global location of the left end of the beam
                  
directionData='C:\Users\luizv\OneDrive\DynamoRevit\Dynamo_visualization\';
ExportResultsBeam(directionData,dim_beams_collection,coordEndBeams,...
                disposition_rebar_beams3sec,nbar_beams_collection,...
                collection_beams_diamIndex,cols_sym_asym_isr,...
                shear_beam_design_collec);
