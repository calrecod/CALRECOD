function [d,qprom]=shearFootings(be,le,qprom,dimCol,pu,d,fc,typeFooting)

%------------------------------------------------------------------------
% Syntax:
% [d,qmax]=shearFootings(be,le,qprom,dimCol,pu,d,fc,typeFooting)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the demand of shear stresses and resistant ones of an
% isolated footing subject to biaxial eccentric actions, considering two 
% mechanisms (punching and beam shearing).
% 
% OUTPUT: d:               effective modified height dimension based on the
%                          critical acting shear demand over the isolated 
%                          footing
%
% INPUT:  qprom:           is the average contact pressure considering the
%                          positive pressures at the four footing's 
%                          corners
%
%         dimCol:          are the cross-section dimensions of the column
%                          that the footing supports (cm)
%
%         be,le:           are the transversal dimensions of the isolated 
%                          footing on plan view (cm)
%
%         pu:              is the axial load reaction from the column (Kg)
%
%         d:               is the effective footing height (cm)
%
%         fc:              is the f'c used for the footing (Kg/cm2)
%
%         typeFooting:     1 - Standard Isolated footing
%                          2 - Bordering Isolated footing
%                          3 - Corner Isolated footing
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

% -----------------------------------------------------------------------
% Shear revision 
% -----------------------------------------------------------------------

%% Punching shear:
vneta=abs(pu)-qprom*dimCol(1)*dimCol(2);
if typeFooting==1
    ashear=d*((dimCol(1)+d)*2+(dimCol(2)+d)*2);
elseif typeFooting==2
    ashear=d*((dimCol(1)+d)+(dimCol(2)+d)*2);
elseif typeFooting==3
    ashear=d*((dimCol(1)+d)+(dimCol(2)+d));
end


vup=vneta/ashear; % punching shear demand
vcp=0.85*fc^0.5; % resistant punching shear

dpunz=d;
if vup>vcp
    while vup>vcp
        dpunz=dpunz+5; % minimum required
        if typeFooting==1
            ashear=dpunz*((dimCol(1)+dpunz)*2+(dimCol(2)+dpunz)*2);
        elseif typeFooting==2
            ashear=dpunz*((dimCol(1)+dpunz)+(dimCol(2)+dpunz)*2);
        end
        vup=vneta/ashear;

    end
elseif vup<vcp
    while vup<vcp
        if dpunz>15
            dpunz=dpunz-5; % minimum required
        end
        if typeFooting==1
            ashear=dpunz*((dimCol(1)+dpunz)*2+(dimCol(2)+dpunz)*2);
        elseif typeFooting==2
            ashear=dpunz*((dimCol(1)+dpunz)+(dimCol(2)+dpunz)*2);
        end
        vup=vneta/ashear;
        if dpunz<=15
            break;
        end
    end
end

%% Revision of shear by flexure (as a beam):
vcrflex=0.5*0.85*fc^0.5; % resistant shear as a beam

% Along the dimension Be
if typeFooting==1
    bp=(be-dimCol(2)-d)*0.5;
elseif typeFooting==2 || typeFooting==3
    bp=(be-dimCol(2)-d);
end
fvi=qprom*bp*le; % Total shear design force

dbe=fvi/(le*vcrflex); % minimum height required as a beam

% Along the dimension Le
if typeFooting==1 || typeFooting==2
    lp=(le-dimCol(1)-d)*0.5;
elseif typeFooting==3
    lp=(le-dimCol(1)-d);
end
fvi=qprom*lp*be; % Total shear design force

dle=fvi/(be*vcrflex); % minimum height required as a beam

%% The max required height considering both mechanisms
vector_d=[dpunz,dbe,dle];
d=max(vector_d);

if mod(d,5)~=0
    d=d+5-mod(d,5); % effectve modified height rounded to a multiply of 5
end

