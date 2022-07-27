function [PC]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)

%------------------------------------------------------------------------
% Syntax:
% [PC]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
%                           rebarAvailable)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the location of the Plastic Center for an asymmetrically
% reinforced concrete cross-section in the axis of reference.
% 
% OUTPUT: PC:                   is the location of the Plastic Center from
%                               the outer most cross-section fiber in the
%                               axis of reference (cm)
%
% INPUT:  b,h:                  cross-section initial dimensions (cm)
%
%         RebarAvailable:       rebar database consisting of an array of 
%                               size [7,3] by default in format: 
%                               [#rebar,diam,unit-weight]
%
%         dispositionRebar:     are the local rebar coordinates
%
%         fy:                   reinforing steel's yield stress (Kg/cm2)
%
%         fdpc:                 is the f'c of concrete reduced by 0.85,
%                               according to the ACI 318 design code
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-07-19
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
                    
    % Plastic Center ____________________________________
    nv=length(rebarTypeslist);
    fd=0;
    f=0;
    for j=1:nv
       av=(rebarAvailable(rebarTypeslist(j),2)^2*pi/4);
       yp=dispositionRebar(j,2);
       d=h/2-yp;
       fd=fd+av*fy*d;
       f=f+av*fy;
    end
    PC=(fdpc*b*(h^2)/2+fd)/(fdpc*b*h+f);