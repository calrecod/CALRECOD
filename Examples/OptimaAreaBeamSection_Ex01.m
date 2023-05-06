% OptimaAreaBeamSection_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the minimum required reinforcement area in tension for a
%    rectangular beam cross-section, according the uniaxial pure flexure 
%    loads
%
%    Note: function SGD1tBeamsISR is the only one required to
%          determine such optimal reinforcement area
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=20; % cross-section width
h=40; % cross-section height
b_rec=4; % concrete cover in the horizontal direction
h_rec=3; % concrete cover in the vertical direction

%% Materials
fc=280; % concrete's compressive strength
E=2e6; % Modulus of Elasticity of the reinforcing steel

%% Additional parameters
duct=3;
factor_fc=0.85; % to compute the reduced f'c

%% Loads
load_conditions=[1 15e5]; % [n-load, Mu]

%% Optimization process
[cbest,bestMr,bestef,best_Area,tbest,h]=SGD1tBeamsISR(b,h,duct,...
            b_rec,h_rec,fc,load_conditions,factor_fc,E,1)
                    
