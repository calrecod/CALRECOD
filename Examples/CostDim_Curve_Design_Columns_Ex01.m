% CostDim_Curve_Design_Columns_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To plot the dimension-cost curves of a reinforced rectangular column
%    cross-section subject to biaxial bending-compression so that the
%    optimal cross-section dimensions can be determined for a given set of
%    load conditions. The function considers the cost of concrete and
%    reinforcing steel.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-26
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-04-16
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% GEOMETRY
height=400; %cm
rec=[4 4]; % [coverx covery] (cm)

%% MATERIALS
fc=280; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
wac=7.8e-3; % unit volume weight of the reinforcing steel (Kg/cm3)

%% LOAD CONDITIONS
load_conditions=[1 15e3 32e5 23e5]; % [nload, Pu, Mx, My] (Kg-cm)

%% Plot options 
plots=1; % for reinforced cross-section plotts (0-No,1-Yes)
graphConvergencePlot=0; % for optimal ISR area convergence (0-No,1-Yes)

%% Additional data
pu_cols=[28.93];

% Ductility demand
ductility=3;

             %f'c % Bombed %Direct Shot (per unit volume)
cost_concrete=[100 2281.22e-6 2266.98e-6;
                150 2401.22e-6 2390.88e-6;
                200 2532.14e-6 2525.28e-6;
                250 2777.96e-6 2845.00e-6;
                300 2939.12e-6 3010.90e-6;
                350 3111.50e-6 3188.50e-6;
                400 3298.16e-6 3380.50e-6;
                450 3499.10e-6 3587.35e-6];

% Unitary cost of concrete per unit volume
nconcrete=length(cost_concrete(:,1));
for i=1:nconcrete-1
    if cost_concrete(i,1)==fc
        unit_cost_conc_cols=cost_concrete(i,2);
        break;
    elseif cost_concrete(i,1)<fc && cost_concrete(i+1,1)>fc
        unit_cost_conc_cols=0.5*(cost_concrete(i,2)+cost_concrete(i+1,2));
        break;
    end
end

%% Cost-dimension curves
[collectionDimCols,collectionISRareaCols,collectionISRcols,...
collectionEffCols]=CostDimCurveOptimDesignCols(height,wac,fc,rec,...
ductility,pu_cols,unit_cost_conc_cols,load_conditions,fy);
%