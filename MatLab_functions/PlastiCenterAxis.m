function [PC]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)

%------------------------------------------------------------------------
% Syntax:
% [PC]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
%                           rebarAvailable)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the location of the Plastic Center for an asymmetrically
% reinforced concrete cross-section in the axis of reference.
% 
% OUTPUT: PC:                   is the location of the Plastic Center from
%                               the outer most cross-section fiber in the
%                               axis of reference
%
% INPUT:  b,h:                  cross-section dimensions
%
%         RebarAvailable:       rebar database consisting of an array of 
%                               size nbars x 2, in format: 
%                               [#rebar,diam]
%
%         dispositionRebar:     are the local rebar coordinates
%
%         fy:                   reinforing steel's yield stress
%
%         fdpc:                 is the f'c of concrete reduced by 0.85,
%                               according to the ACI 318 design code
%
%         rebarTypeslist:       is the list containing the rebar diameter 
%                               index corresponding to its location in the
%                               array RebarAvailable, for each cross-
%                               section's edge
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
                    
% Plastic Center ________________________________________________________
nv=length(rebarTypeslist);
fd=0;
f=0;
for j=1:nv
   av=(rebarAvailable(rebarTypeslist(j),2))^2*pi/4;
   yp=dispositionRebar(j,2);
   d=h/2-yp;
   fd=fd+av*fy*d;
   f=f+av*fy;
end
PC=(fdpc*b*(h^2)/2+fd)/(fdpc*b*h+f);