% Design_Footings_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally the reinforcement of an isolated footing's both
%    transversal cross-section 
%
%    Note: function designDimFootings is used to dimension the transversal
%          cross section dimensions of the footing, based on the bearing
%          load capacity of the soil and the load pressures
%
%          function isrFootings is the one required to optimally design the
%          reinforcement. Two options are available for reinforcement:
%          either a pure ISR or with rebar
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clear all
clc

%% Geometry data
rec=5; % concrete cover in all directions
hefoot=25; % initial footing height dimensions (cm)
d=hefoot-rec; % effective footing's height (cm)
bc=30; hc=30; % supporting column's dimensions
dimCol=[bc,hc]; 

%% Materials
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
fc=250; % concrete's compressive strength f'c (Kg/cm2)
wac=7.8e-3; % unit volume weight of the reinforcing steel (Kg/cm3)

%% Additional data
qadm=1.5; % Admissible bearing load of soil (Kg/cm2)
FS=1.3; % safety factor for the soil's bearing load capacity 
qu=qadm/FS;
pu_steel_footings=26.75; % unit construction cost of reinforcement assembly
                         % (Per kilogram of reinforcing steel)             

typeFoot=3; % Type of footing (1, 2 or 3) - see documentation
sepMinRebar=15; % Minimum rebar separation considered (cm) greater or 
                % equal than the min absoulute give by code: 4/3*0.75 in

cols_sym_asym_isr="Rebar"; % ''ISR'' or anything else (''Symmetric'', ''Rebar'')
ductility=3; % to set the MAX reinforcement quantity (1, 2 or 3)

optimConv=1; % to visualize the optimization convergence plot
PlotRebarDesign=1; % To visualize the optimal rebar design layouts

%% Loads
load_cols=[1 -14.59e3 -1.2e5 1.13e5]; % [n-load, Pu, Mux, Muy] Kg,cm
pu=load_cols(1,2);

%% Rebar data
% Commercially available rebar
                %type diam 
RebarAvailable=[4 4/8.*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
%% Design of footing's transversal dimensions
% Note: in case it is required to design transversal dimensions
[be,le,contact_pressure]=designDimFootings(pu,qu,dimCol,hefoot,rec,...
                        typeFoot)

dim_zap=[be le];

%% Optimization design process
[hmodif,m_max_eje,barDispositionFootings,arrangement_bar_footings,...
nbars_footings,AcBar,bestCost_elem,list_ef_footings,list_mr_footings]=...
isrFootings(pu_steel_footings,hefoot,dim_zap(1),dim_zap(2),...
rec,fc,fy,load_cols,dimCol,RebarAvailable,cols_sym_asym_isr,...
ductility,optimConv,PlotRebarDesign,typeFoot,sepMinRebar,wac);

%% Exporting design results for visualization
dimensionFootingCollection=[be,le,hmodif,rec]
coordBaseFooting=[0,0,0]; % global position of the footing

directionData=[];
ExportResultsIsolFootings(directionData,barDispositionFootings,...
    dimensionFootingCollection,nbars_footings,arrangement_bar_footings,...
    coordBaseFooting,cols_sym_asym_isr)
