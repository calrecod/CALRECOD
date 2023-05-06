function [disposition_rebar,separation_hor1,separation_hor2,...
    separation_ver1,separation_ver2]=dispositionRebarAsymmetric(b,...
    h,rec,nv,number_rebars_sup,number_rebars_inf,number_rebars_left,...
    number_rebars_right,RebarAvailable,op1,op2,op3,op4)

%------------------------------------------------------------------------
% Syntax:
% [disposition_rebar,separation_hor1,separation_hor2,...
%   separation_ver1,separation_ver2]=dispositionRebarAsymmetric(b,...
%   h,rec,nv,number_rebars_sup,number_rebars_inf,number_rebars_left,...
%   number_rebars_right,RebarAvailable,op1,op2,op3,op4)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local coordinates of an asymmetric rebar option.
% 
% OUTPUT: disposition_rebar:    are the local coordinates of the optimal rebar option
%
%         separation_hor1,
%         separation_hor2,
%         separation_ver1,
%         separation_ver2:      resultant rebar separation to be compared 
%                               with the minimum one (upper, lower, left
%                               right boundary), respectively
%
% INPUT:  b,h:                  cross-section initial dimensions
%
%         rec:                  are the concrete cover values for each axis
%                               direction of the cross-section
%
%         RebarAvailable:       rebar database consisting of an array of 
%                               size: n# x 3, by default in format: 
%                               [#rebar,diam,unit-weight]
%
%         number_rebars_sup,
%         number_rebars_inf,
%         number_rebars_left,
%         number_rebars_right:  number of rebars to placed on each of the 
%                               cross-section boundaries
%
%         op1,op2,op3,op4:      resultant types of rebar for each of the 
%                               four cross-section boundaries (upper 
%                               boundary, lower boundary, left side and 
%                               right side, respectively)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bprima=b-2*rec(1);
hprima=h-2*rec(2);

dv1=RebarAvailable(op1,2);
dv2=RebarAvailable(op2,2);
dv3=RebarAvailable(op3,2);
dv4=RebarAvailable(op4,2);

ndiam=length(RebarAvailable(:,1));
disposition_rebar=zeros(nv,2);
dv=RebarAvailable(ndiam,2);
if number_rebars_sup>=2
    
    separation_hor1=round((bprima-((number_rebars_sup)*dv1))/...
        (number_rebars_sup-1),1);
    
    separation_ver2=round(((hprima)-((number_rebars_right+2)*dv4))/...
    (number_rebars_right+1),1);

    % Rebar disposition-----------------------------------------------
    % rebars disposition horizontally (superior)
    for j=1:number_rebars_sup

        % rebar coordinates in the superior part
        disposition_rebar(j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
            (separation_hor1+dv1);
        disposition_rebar(j,2)=hprima*0.5-(dv*0.25);
    end
    
    % rebar disposition in the vertical right part
    for j=1:number_rebars_right

        % right side
        disposition_rebar(j+number_rebars_sup+number_rebars_inf+...
            number_rebars_left,2)=hprima*0.5-(dv4)*0.5-(separation_ver2+dv4*0.5)-(j-1)*(separation_ver2+dv4);

        disposition_rebar(j+number_rebars_sup+number_rebars_inf+...
            number_rebars_left,1)=bprima*0.5-(dv*0.5);
    end

else
    separation_ver2=round((hprima-((number_rebars_right+1)*dv4))/...
    (number_rebars_right),1);

    % REBAR DISPOSITION-----------------------------------------------
    % rebar disposition horizontally (superior) only one rebar
    for j=1:number_rebars_sup

        % rebar coordinates in the superior part
        disposition_rebar(j,1)=-bprima*0.5+(dv)*0.5;
        disposition_rebar(j,2)=hprima*0.5-(dv*0.5);
    end
    
    % rebar disposition over right side
    for j=1:number_rebars_right
        % right side
        disposition_rebar(j+number_rebars_sup+number_rebars_inf+...
            number_rebars_left,2)=hprima*0.5-(j-1)*(separation_ver2+dv4);

        disposition_rebar(j+number_rebars_sup+number_rebars_inf+...
            number_rebars_left,1)=bprima*0.5-(dv*0.5);
    end
end

if number_rebars_inf>=2
    separation_hor2=round((bprima-((number_rebars_inf)*dv2))/...
        (number_rebars_inf-1),1);
    
    separation_ver1=round((hprima-((number_rebars_left+2)*dv3))/...
    (number_rebars_left+1),1);
    
    % REBAR DISPOSITION-----------------------------------------------
    % rebar coordinates on the inferior part
    for j=1:number_rebars_inf
        
        disposition_rebar(number_rebars_sup+j,1)=-bprima*0.5+...
            (dv)*0.5+(j-1)*(separation_hor2+dv2);
        disposition_rebar(number_rebars_sup+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % left side
    for j=1:number_rebars_left

        disposition_rebar(j+number_rebars_sup+number_rebars_inf,2)=...
            hprima*0.5-(dv3)*0.5-(separation_ver1+dv3*0.5)-(j-1)*(separation_ver1+dv3);

        disposition_rebar(j+number_rebars_sup+number_rebars_inf,1)=...
            -bprima*0.5+(dv*0.5);
    end
else
    separation_ver1=round((hprima-((number_rebars_left+1)*dv3))/...
    (number_rebars_left),1);

    % REBAR DISPOSITION-----------------------------------------------
    % rebar coordinates over the inferior part
    for j=1:number_rebars_inf
        
        disposition_rebar(number_rebars_sup+j,1)=bprima*0.5-...
            (dv)*0.5;
        disposition_rebar(number_rebars_sup+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % left side
    for j=1:number_rebars_left

        disposition_rebar(j+number_rebars_sup+number_rebars_inf,2)=...
            -hprima*0.5+(dv3)*0.5+(j-1)*(separation_ver1+dv3);

        disposition_rebar(j+number_rebars_sup+number_rebars_inf,1)=...
            -bprima*0.5+(dv*0.5);

    end

end