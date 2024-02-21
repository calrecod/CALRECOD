% CostDim_Curve_Design_Beams_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To plot the dimension-cost curves of a reinforced rectangular beam
%    cross-section subject to pure flexure so that the
%    optimal cross-section dimensions can be determined for a given set of
%    load conditions. The function considers the cost of concrete and
%    reinforcing steel.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-26
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% Geometry
span=500; % beam's length (cm)
                
h_rec=[4]; % vertical concrete cover  
b_rec=3; % lateral concrete cover

%% Material
fc=280; % Concrete compressive strength Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
wac=7.8e-3; % unit volume weight of the reinforcing steel (Kg/cm3)

%% Load conditions
load_conditions=[1 -23.0e5 16.0e5 -20.0e5]; % Kg-cm

%% Additional data
pu_beams_steel=38.6; % unit assembly cost of steel reinforcement
duct=1; % high ductility demand

             %f'c % Bombed %Direct Shot (per unit volume)
cost_concrete=[100 2281.22e-6;
                150 2401.22e-6;
                200 2532.14e-6;
                250 2777.96e-6;
                300 2939.12e-6;
                350 3111.50e-6;
                400 3298.16e-6;
                450 3499.10e-6];

%% Unit cost of concrete per unit volume
nconcrete=length(cost_concrete(:,1));
for i=1:nconcrete-1 % a loop for seraching the corresponding unit cost
                    % according to the f'c
    if cost_concrete(i,1)==fc
        unit_cost_conc_beams=cost_concrete(i,2);
        break;
    elseif cost_concrete(i,1)<fc && cost_concrete(i+1,1)>fc
        unit_cost_conc_beams=0.5*(cost_concrete(i,2)+cost_concrete(i+1,2));
        break;
    end
end

%% Cost-Dim curve design
[collectionDimBeams,collectionISRareaBeams,collectionISRbeams,...
    collectionEffBeams]=CostDimCurveOptimDesignBeam(span,wac,fc,b_rec,...
    h_rec,duct,pu_beams_steel,unit_cost_conc_beams,load_conditions,fy);
%