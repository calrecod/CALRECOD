% StrucEfficiency_TBeamsRebar_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To analise the structural efficiency of a rebar distribution over a
%    beam element of T cross-section 
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-04-09
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

%% Geometry

bp=20; % web width (cm) 
ht=30; % total height (cm)
ba=60; % flange width (cm) 
ha=12; % flange height or thickness (cm)
span=500; % cm

cover=3; % lateral concrete cover

%% Materials
fc=250; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)

fdpc=fc*0.85; % reduced f'c
beta1=0.85;

%% Load conditions
load_conditions=[1 700000.0];

%% Rebar data
% Database of the commercially available rebar
rebarAvailable=[3 3/8*2.54;
                4 4/8.*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
% Rebar coordinates
dispositionRebar=[-7 -12;
                   0 -12;
                   7 -12];
rebarType=[4;4;4];
            
As=sum(rebarAvailable(rebarType,2).^2*pi./4);

%% Analysis of efficiency
[Eff,Mr,cx]=EfRebarTBeams(load_conditions,bp,ht,ba,ha,span,fdpc,...
                           rebarType,rebarAvailable,cover,beta1,...
                           dispositionRebar)
                       
