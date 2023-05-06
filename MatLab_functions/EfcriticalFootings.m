function [maxef,mr]=EfcriticalFootings(bz,fdpc,actension,d,fy,mu_real_axis)

%------------------------------------------------------------------------
% Syntax:
% [maxef,mr]=EfcriticalFootings(bz,fdpc,actension,d,fy,mu_real_axis)
%
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
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
%         d:                    effective footing height h-cover
%
%         actension:            is the quantity of rebar area in tension
%
%         fdpc:                 is the f'c used reduced by the 0.85 factor 
%                               defined by code
%
%         fy:                   is the yield strength of the rebar steel
%
%         mu_real_axis:         is the demanding bending moment
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

c=actension*fy/(0.85*bz*fdpc);
if c>(d/(0.004/0.003+1))
    c=d/(0.004/0.003+1);
end
    a=0.85*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=actension*fy;
mr=0.9*(T*(d-a/2));

maxef=mu_real_axis/mr;
