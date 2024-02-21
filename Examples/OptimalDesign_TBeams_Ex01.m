% OptimalDesign_TBeams_Ex01
%------------------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to savings in reinforcing area)
%    a beam element of T cross-section 
%
%------------------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
clc
clear all

%% Geometry
bp=20; % web width (cm) 
ht=30; % total height (cm)
ba=60; % flange width (cm) 
ha=12; % flange height or thickness (cm)
span=500; % beam's length (cm)

cover=3; % concrete cover
d=ht-cover; % effective cross-section's height

%% Material
fc=250; % concrete's compressive strength Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
Es=2.0e6;
ffc=0.85;
fdpc=fc*ffc;

%% Load conditions
load_conditions=[1 6.55e5]; % [n-load, Mu]

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
            
wac=7.8e-3; % unit volume weight of the reinforcing steel (Kg/cm3)

%% Additional data:
duct=3;

cols_sym_asym_isr="Rebar"; % ''Rebar'' or ''ISR''

% Plot options:
rebarDesignPlots=1; % for reinforced cross-section plots (0-No,1-Yes)
graphConvergencePlot=1; % to visualize the ISR optimization convergence

puTbeams=41.6; % unit construction assembly cost of steel reinforcement

%% Optimization through the ISR analogy:
[cbest,bestMr,bestef,best_Area,tbest]=SGD1tTBeamsISR(bp,ht,ba,ha,span,duct,...
      cover,fc,load_conditions,ffc,Es,graphConvergencePlot);
        
% Compute cost with the ISR:
bpp=bp-2*cover;
tmin=(0.7*sqrt(fdpc)/fy*(bp*(d-ha)+ha*ba))/bpp; % min ISR's width in 
                                                % compression
t2Best=[tbest,tmin]; % ISR's widths in tension and compression

%% Rebar design optimization:
if cols_sym_asym_isr~="ISR"
        
    [sepbarsRestric,cbest,bestBarDisposition,bestCostRebar,barTypes1Comp,...
        barTypesTen,ef,bestMr,areaRebar]=ISR1tRebarTBeamsOptim(bp,ht,...
        ba,ha,fc,cover,load_conditions,t2Best,puTbeams,span,rebarAvailable,...
        wac);
end
 
%% Plotting reinforced cross-section:
if rebarDesignPlots==1
    TbeamReinforcedSection(bp,ht,ba,ha,bestBarDisposition,...
        barTypes1Comp,barTypesTen);
end
