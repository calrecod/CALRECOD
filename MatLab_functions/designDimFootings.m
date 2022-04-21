function [be,le,contact_pressure]=designDimFootings(pu,qu,dimCol,...
    hfooting,rec)

%------------------------------------------------------------------------
% Syntax:
% [be,le,contact_pressure]=designDimFootings(pu,qu,dimCol,...
%   hfooting,rec)
%
%------------------------------------------------------------------------
% PURPOSE: To design the transversal dimensions of a rectangular isolated 
% footing, based on the acting vertical reaction from the intersecting 
% column and the admissible load of soil considering the Safety Design Factor,
% as well as the cross-section dimensions of the intersecting column, so 
% that the footing dimension may be at least 40 cm wider than such column
% cross-section dimensions.
% 
% OUTPUT: be,le:                transversal isolated footing dimensions
%
%         contact_pressure:     is the resulting contact pressure from the
%                               soil (less or equal than qu)
%
% INPUT:  pu:                   is the vertical reaction from he supporting 
%                               column onto the footing
%
%         qu                    is equal to FS(q{adm})
%
%         dimCol:               column cross-section dimensions [b,h]
%
%         hfooting:             is the height or width of the isolated
%                               footing
%
%         rec:                  is the concrete cover (cm)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

lebe=sqrt(pu/qu);
be=lebe;
le=lebe;
if (be-dimCol(1,2))*0.5<((hfooting-rec)*0.5+10)
    be=dimCol(1,2)+20+(hfooting-rec);
end

if (le-dimCol(1,1))*0.5<((hfooting-rec)*0.5+10)
    le=dimCol(1,1)+20+(hfooting-rec);
end
contact_pressure=pu/(be*le);
