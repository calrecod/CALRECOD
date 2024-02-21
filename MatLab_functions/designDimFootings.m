function [be,le,contactPressure]=designDimFootings(pu,qmax,dimCol,...
    hfooting,rec,typeFooting)

%------------------------------------------------------------------------
% Syntax:
% [be,le,contactPressure]=designDimFootings(pu,qadm,dimCol,...
%   hfooting,rec,typeFooting)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To design the transversal dimensions of a rectangular isolated 
% footing, based on the acting vertical reaction from the intersecting 
% column and the admissible load of soil by considering the cross-section 
% dimensions of the intersecting column, so that the footing dimension 
% may be at least 40 cm wider than such column cross-section dimensions.
% 
% OUTPUT: be,le:                transversal isolated footing dimensions
%
%         contactPressure:      is the resulting contact pressure from the
%                               soil (less or equal than qu)
%
% INPUT:  pu:                   is the vertical reaction from the supporting 
%                               column onto the footing
%
%         qmax                  The max bearing soil capacity
%
%         dimCol:               column cross-section dimensions [b,h]
%
%         hfooting:             is the height or thickness of the isolated
%                               footing
%
%         rec:                  is the concrete cover
%   
%         typeFooting:          1 - Standard Isolated footing
%                               2 - Bordering Isolated footing
%                               3 - Corner Isolated footing
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

% Dimensions of the column ----------------------------------------------
bc=dimCol(1);
hc=dimCol(2);
% -----------------------------------------------------------------------

lebe=sqrt(abs(pu)/qmax); % Considering initially a square footing
be=lebe;
le=lebe;

if typeFooting==1 % Standard Isolated Rectangular Footing
    if (be-hc)*0.5<((hfooting-rec)*0.5+20)
        be=hc+40+(hfooting-rec);
    end

    if (le-bc)*0.5<((hfooting-rec)*0.5+20)
        le=bc+40+(hfooting-rec);
    end
elseif typeFooting==2 % Bordering Isolated Rectangular Footing
    if (be-hc)<((hfooting-rec)*0.5+40)
        be=hc+40+(hfooting-rec)*0.5;
    end

    if (le-bc)*0.5<((hfooting-rec)*0.5+20)
        le=bc+40+(hfooting-rec)*0.5;
    end
elseif typeFooting==3 % Corner Isolated Rectangular Footing
    if (be-hc)<((hfooting-rec)*0.5+40)
        be=hc+40+(hfooting-rec)*0.5;
    end

    if (le-bc)<((hfooting-rec)*0.5+40)
        le=bc+40+(hfooting-rec)*0.5;
    end
end
be=be-mod(be,5)+5; % The dimensions are rounded to a multply of 5
le=le-mod(le,5)+5;

contactPressure=abs(pu)/(be*le);
