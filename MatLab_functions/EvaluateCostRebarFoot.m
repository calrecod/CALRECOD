function cost = EvaluateCostRebarFoot(be,RebarAvailable,...
                                arrangement1,arrangement2,pu)

%------------------------------------------------------------------------
% Syntax:
% cost = EvaluateCostRebarFoot(be,RebarAvailable,...
%                              arrangement1,arrangement2,pu)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the construction cost of a rebar option for an 
% isolated footing transversal cross-section.
% 
% OUTPUT: cost:                 is the construction cost of a rebar option
%                               over an isolated footing transversal
%                               cross-section
%
% INPUT:  be:                   longitudinal dimension perpendicular to the 
%                               transversal cross-section in analysis
%                               (length of rebars)
%
%         arrangement1,
%         arrangement2:         are the rebar types used for the steel in
%                               tension and compression, respectively
%
%         pu:                   is the unit construction assembly cost of 
%                               rebars in isolated footings: units $/Kg. 
%                               The unit cost is considered using an average
%                               value of all rebar types's assembly 
%                               performances (assuming that in an isolated
%                               footing as much as four different types of
%                               rebar may be placed)
%
%         RebarAvailable:       is the array containing the available 
%                               commercial rebar database: size = [7,5]
%                               by default
%                   _________________________________________
%                   [noption,# rebar,diam,area,lineal-weight]
%                   -----------------------------------------
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nvHor1=length(arrangement1);
nvHor2=length(arrangement2);

cost=(nvHor1*RebarAvailable(arrangement1(1),4)*0.0001*7800*be*0.01+...
       nvHor2*RebarAvailable(arrangement2(1),4)*0.0001*7800*be*0.01)*pu;
