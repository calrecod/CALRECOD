% Design_Beams_Shear_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design the transversal rebar along a rectangular beam
%
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------
clc
clear all

%% Geometry data
span=500; % cm
b=30; % width (cm)
h=60; % height (cm)

rec=3; % lateral concrete cover
rho=0.0029; % average steel rebar percentage over the beam cross-sections

%% Material data
fc=280; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)

%% Shear forces along the length of the beam
np=7; % number of points for the distribution of the shear forces
shear_beam=linspace(12000,-22000,np); % This is the vector containing the distributed
									  % shear forces along the beam span

[s1,s2,s3,d1,d3]=shearDesignBeams(span,b,h,rec,rho,fc,fy,...
                          shear_beam)
