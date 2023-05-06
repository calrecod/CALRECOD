% Efficiency_Beam_RebarDesign_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the structural efficiency of a
%    reinforced rectangular beam cross-section.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2023-03-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-04-16
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=25; % cross-section width
h=50; % cross-section height

b_rec=3; % concrete cover on each direction
h_rec=3;

d=h-h_rec; % effective cross-section height

%% Materials
fc=250; % concrete's compressive strength
factor_fc=0.85; % reduction factor for the f'c
fdpc=factor_fc*fc; % reduced f'c
E=2e6; % Modulus of elasticity of the reinforcing steel
fy=4200; % Yield stress of the reinforcing steel

%% Additional parameters
duct=3; % ductility demand level

%% Loads
load_conditions=[1 -5.26e5]; % [n-load, Mu] (Kg-cm)

%% Rebar data
% Commercially available rebar diameters
                %type diam  
rebarAvailable=[4 4/8.*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
% Distribution of rebars over the cross-section
dispositionRebar=[-9 21;
                    9 21;
                    -9 -21;
                    0 -21;
                    9 -21];
                
RebarIndexTen=[1;1;1]; % rebar diameters to use for the reinforcement
                    % in tension (indices from the "rebarAvailable" array)
                    
RebarIndexCom=[1;1]; % rebar diameters to use for the reinforcement
                  % in compression (indices from the "rebarAvailable" array)

%% Additional design information of interest
ast=sum(rebarAvailable(RebarIndexTen,2).^2.*pi./4);

astotal=ast+asc % Total rebar area
amax=0.025*b*d % Max allowed rebar area by code

%% Structural efficiency
[maxef,Mrv,c]=EfcriticalRebarbeams(load_conditions,b,E,fdpc,RebarIndexTen,...
    RebarIndexCom,rebarAvailable,d,h_rec,0.85,dispositionRebar)

