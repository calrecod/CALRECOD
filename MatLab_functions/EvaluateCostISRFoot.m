function cost = EvaluateCostISRFoot(be,rec,act,acmin,pu_steel_footings,wac)

%------------------------------------------------------------------------
% Syntax:
% cost = EvaluateCostISRFoot(be,rec,act,acmin,pu,wac)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the estimated construction cost given an ISR data.
% 
% OUTPUT: cost:                 is the estimated construction cost of an ISR 
%                               over an isolated footing transversal 
%                               cross-section
%
% INPUT:  be:                   longitudinal dimension perpendicular to the
%                               transversal cross-section in analysis
%                               (length of rebars)
%
%         act:                  is the steel area reinforcement in tension
%                               over the transversal cross-section (lower 
%                               boundary)
%
%         acmin:                is the minimum steel area reinforcement by 
%                               temperature over the zone in compression of 
%                               the transversal cross-section (upper 
%                               boundary)
%
%         rec:                  is the concrete cover
%
%         pu:                   is the unit construction assembly cost of
%                               rebars in isolated footings: units $/Kg. 
%                               The unit cost is considered using an average
%                               value of all rebar types's assembly 
%                               performances (assuming that in an isolated
%                               footing as much as four different types of 
%                               rebar may be placed)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

cost=act*(be-2*rec)*wac*pu_steel_footings+...
     acmin*(be-2*rec)*wac*pu_steel_footings;
                    