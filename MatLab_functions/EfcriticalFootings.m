function [maxef,mr]=EfcriticalFootings(bz,fdpc,actension,d,fy,mu_real_axis)

%------------------------------------------------------------------------
% Syntax:
% [maxef,mr]=EfcriticalFootings(bz,fdpc,actension,d,fy,mu_real_axis)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the structural efficiency of a rebar option, 
% considering only the steel in tension for a transversal cross-section of
% an isolated footing in which the steel in compression is designed only
% by temperature.
% 
% OUTPUT: maxef,mr:             is the structural efficiency of the designed
%                               transversal cross-section and the resisting
%                               bending moment, respectively
%
% INPUT:  bz:                   width dimension of the transversal cross
%                               section in analysis
%
%         d:                    effective footing height h-cover (cm)
%
%         actension:            is the quantity of rebar area in tension
%
%         fdpc:                 is the f'c used reduced by the 0.85 factor 
%                               defined by code (Kg/cm2)
%
%         fy:                   is the yield strength of the rebar steel
%                               units (Kg/cm2)
%
%         mu_real_axis:         is the demanding bending moment (Kg-cm)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

c=actension*fy/(0.85*bz*fdpc);
    a=0.85*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=actension*fy;
mr=0.9*(T*(d-a/2))*1e-5;

maxef=mu_real_axis*0.00001/mr;
