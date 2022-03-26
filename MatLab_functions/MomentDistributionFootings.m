function [m]=MomentDistributionFootings(sigma1,sigma2,qmax1,qmax2,...
            a1,a2,dimp,dimfoot)

%------------------------------------------------------------------------
% Syntax:
% [m]=MomentDistributionFootings(sigma1,sigma2,qmax1,qmax2,...
%            a1,a2,dimp,dimfoot)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the effective bending moment acting at the transversal
% cross-sections of the footing, given the max unit contact pressures at 
% each of the four footing's boundaries.
% 
% OUTPUT: m:                  effective bending moment acting over a 
%                             transversal footing cross-section
%
% INPUT:  sigma1,sigma2:      are magnitudes of the contact pressures at 
%                             the other end of both effective unit contact
%                             pressures at each side of the axis direction 
%                             in question. (see Documentation)
%
%         qmax1,qmax2:        are the max contact pressures at both sides
%                             of the axis direction of anlaysis. (See
%                             Documentation)
%
%         a1,a2:              are the effective unit contact pressures at
%                             both sides of the axis direction of analysis.
%                             (see Documentation)
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

%---- Para a_13----%
a1x=sigma1*dimp*(dimp*0.5)+0.5*dimp*(qmax1-sigma1)*(2/3*dimp);
x12=a1x/a1;

%---- Para a_13----%
a2x=sigma2*dimp*(dimp*0.5)+0.5*dimp*(qmax2-sigma2)*(2/3*dimp);
x34=a2x/a2;

%----- x prom -----%
xt=(x12+x34)*0.5;
aprom=(a1+a2)*0.5;

fcorte=aprom*dimfoot;

m=fcorte*xt;
