function [BS,CFA]=BuildabilityScoreRebarCols(irebar4,nb4,wunb,wnd)
%------------------------------------------------------------------------
% Syntax:
% [BS,CFA]=BuildabilityScoreRebarCols(irebar4,nb4,wunb,wnd)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: None
%
%------------------------------------------------------------------------
% PURPOSE: To determine the Buildability Score BS and the Complexity Factor
% of Assembly CFA of a rebar design in a rectangular concrete column.
% 
% OUTPUT: BS:              Buildability Score
%
%         CFA:             Complexity Factor of Assembly 
%
% INPUT:  irebar4:         are the rebar indeces of each rebar diameter of
%                          each cross-section's edge from the commercially
%                          available rebar table
%
%         nb4:             is the vector containing the number of rebar on
%                          each cross-section's edge. Size: 1 x 4
%
%         wunb,wnd:        weight factors for the uniformity of number of
%                          rebars (UNB variable) and the Number of
%                          Diameters (ND variable)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-11-11
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

if nb4(1)>nb4(2)
    nb1=nb4(2);
    nb2=nb4(1);
elseif nb4(1)<=nb4(2)
    nb1=nb4(1);
    nb2=nb4(2);
end

if nb4(3)>nb4(4)
    nb3=nb4(4);
    nb44=nb4(3);
elseif nb4(3)<=nb4(4)
    nb3=nb4(3);
    nb44=nb4(4);
end

nb12=nb1/nb2;

if nb44==0
    nb34=1;
else
    nb34=nb3/nb44;
end
unb=0.5*(nb12+nb34);
nd=1;
for i=1:3
    if irebar4(i)~=irebar4(i+1) && nb4(i+1)~=0
        nd=nd+1;
    end
end
BS=unb^wunb+1/nd^wnd;
CFA=BS/2;  