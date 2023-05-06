function [newdispositionRebar,newCoordCorners,newDepthCP]=rotReCol2(Mux,...
                        Muy,dispositionRebar,b,h,CPaxis)
                            
%------------------------------------------------------------------------
% Syntax:
% [newdispositionRebar,newCoordCorners,newDepthCP]=rotReCol2(Mux,...
%  Muy,dispositionRebar,b,h,CPaxis)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
%------------------------------------------------------------------------
% PURPOSE: To rotate a rectangular reinforced concrete column cross-section
% according to a pair of bending moments. The function computes the new
% local coordinates of the rebar in the rotated system of reference, the
% cross-section corners and the Plastic Center location.
% 
% OUTPUT: newdispositionRebar:  are the new local rebar coordinates over
%                               the rotated rectangular cross-section 
%
%         newCoordCorners:      are the new local coordinates of the
%                               rectangular cross-section corners in the
%                               rotated system of reference
%
%         newDpethCP:           is the new depth of the reinforced 
%                               cross-section Plastic Center with respect
%                               to the axis X' and Y' (see doc.)
%
% INPUT:  Mux,Muy:              is the pair of bending moments in the
%                               original non-rotated system of reference
%
%         dispositionRebar:     are the original local rebar coordinates
%                               over the non-rotated rectangular 
%                               cross-section
%
%         b,h:                  are the cross-section dimensions (width and
%                               height, respectively)
%
%         CPaxis:               are the original depths of the
%                               cross-section Plastic Center
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

alpha=rad2deg(atan(Mux/Muy));
if Muy<=0 && Mux>=0 || Muy<=0 && Mux<=0
    gamma=(90-alpha)+180;
elseif Muy>=0 && Mux<=0 || Muy>=0 && Mux>=0
    gamma=90-alpha; % This is the angle for the section to be rotated at
                    % so that the resultant moment of Mx and My is now Mx'
end
nb=length(dispositionRebar(:,1));

%% Rotation of the cross-section corners
coordOrigSec=[0.5*b 0.5*h;
             -0.5*b 0.5*h;
             -0.5*b -0.5*h;
              0.5*b -0.5*h];

% New coordinates of the cross-section corners
newCoordCorners(:,1)=cos(deg2rad(gamma))*coordOrigSec(:,1)+...
                     sin(deg2rad(gamma))*coordOrigSec(:,2);
newCoordCorners(:,2)=-sin(deg2rad(gamma))*coordOrigSec(:,1)+...
                      cos(deg2rad(gamma))*coordOrigSec(:,2);

%% Rotation of the rebars

for i=1:nb
    % New coordinates of each rebar
    newdispositionRebar(i,1)=dispositionRebar(i,1)*cos(deg2rad(gamma))+...
                             dispositionRebar(i,2)*sin(deg2rad(gamma));
    newdispositionRebar(i,2)=-dispositionRebar(i,1)*sin(deg2rad(gamma))+...
                             dispositionRebar(i,2)*cos(deg2rad(gamma));
end

%% Rotation of the Plastic Center location

OrCPaxisCoord=[0.5*h-CPaxis(1),0.5*b-CPaxis(2)];

% New coordinates of the cross-section CP
newDepthCP(:,1)=max(newCoordCorners(:,2))-...
                (OrCPaxisCoord(1,1)*cos(deg2rad(gamma))+...
                OrCPaxisCoord(1,2)*sin(deg2rad(gamma)));
            
newDepthCP(:,2)=min(newCoordCorners(:,1))-...
                (-OrCPaxisCoord(1,1)*sin(deg2rad(gamma))+...
                OrCPaxisCoord(1,2)*cos(deg2rad(gamma)));
            
            