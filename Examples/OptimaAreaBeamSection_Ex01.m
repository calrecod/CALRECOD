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

b=20;
h=40;
duct=3;
b_rec=4;
h_rec=3;
fc=280;
load_conditions=[1 15];
factor_fc=0.85;
E=2e6;
plotOptimConv=1;
[cbest,bestMr,bestef03,best_Area03,tbest,h]=SGD1tBeamsISR(b,h,duct,...
            b_rec,h_rec,fc,load_conditions,factor_fc,E,plotOptimConv)
                    
