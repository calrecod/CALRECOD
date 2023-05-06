function [m]=MomentDistributionFootings(qmax,dimp,dimfoot)

%------------------------------------------------------------------------
% Syntax:
% [m]=MomentDistributionFootings(qmax,dimp,dimfoot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the effective bending moment acting at the transversal
% cross-sections of the footing, given the max unit contact pressures at 
% each of the four footing's boundaries.
% 
% OUTPUT: m:                  effective bending moment acting over a 
%                             transversal footing cross-section
%
% INPUT:  qmax:               is the max contact pressure. (See
%                             Documentation)
%
%         dimp:               is the effective dimension over which the 
%                             flexure stress acts
% 
%         dimfoot:            is the footing dimension over the axis 
%                             direction in analysis
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fshear=qmax*dimp*dimfoot;
xt=0.5*dimp;
m=fshear*xt;
