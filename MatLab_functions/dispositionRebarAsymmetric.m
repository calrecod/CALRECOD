function [disposition_rebar,separation_hor1,separation_hor2,...
    separation_ver1,separation_ver2]=dispositionRebarAsymmetric(b,...
    h,sepMin,rec,nv,number_rebars_sup,number_rebars_inf,number_rebars_left,...
    number_rebars_right,RebarAvailable,op1,op2,op3,op4)

%------------------------------------------------------------------------
% Syntax:
% [disposition_rebar,separation_hor1,separation_hor2,...
%   separation_ver1,separation_ver2]=dispositionRebarAsymmetric(b,...
%   h,sepMin,rec,nv,number_rebars_sup,number_rebars_inf,number_rebars_left,...
%   number_rebars_right,RebarAvailable,op1,op2,op3,op4)
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
%         sepMin:               is the minimum rebar separation constriction 
%                               (see Documentation)
%
%         rec:                  are the concrete cover values for each axis
%                               direction of the cross-section
%
%         RebarAvailable:       rebar database consisting of an array of 
%                               size [7,3] by default in format: 
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

disposition_rebar=zeros(nv,2);
dv=RebarAvailable(7,2);
if number_rebars_sup(op1)>=2
    
    separation_hor1=round((bprima-((number_rebars_sup(op1))*dv1))/...
        (number_rebars_sup(op1)-1),1);
    
    separation_ver2=round(((hprima)-((number_rebars_right(op4)+2)*dv4))/...
    (number_rebars_right(op4)+1),1);

    % Rebar disposition________________________________________________
    %__________________________________________________________________
    % rebars disposition horizontally (superior)
    for j=1:number_rebars_sup(op1)

        % rebar coordinates in the superior part
        disposition_rebar(j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
            (separation_hor1+dv1);
        disposition_rebar(j,2)=hprima*0.5-(dv*0.25);
    end
    
    % rebar disposition in the vertical right part
    for j=1:number_rebars_right(op4)

        % right side
        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2)+...
            number_rebars_left(op3),2)=hprima*0.5-(dv4)*0.5-(separation_ver2+dv4*0.5)-(j-1)*(separation_ver2+dv4);

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2)+...
            number_rebars_left(op3),1)=bprima*0.5-(dv*0.5);
    end

else
    separation_ver2=round((hprima-((number_rebars_right(op4)+1)*dv4))/...
    (number_rebars_right(op4)),1);

    % REBAR DISPOSITION________________________________________________
    %__________________________________________________________________
    % rebar disposition horizontally (superior) only one rebar
    for j=1:number_rebars_sup(op1)

        % rebar coordinates in the superior part
        disposition_rebar(j,1)=-bprima*0.5+(dv)*0.5;
        disposition_rebar(j,2)=hprima*0.5-(dv*0.5);
    end
    
    % rebar disposition over right side
    for j=1:number_rebars_right(op4)
        % right side
        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2)+...
            number_rebars_left(op3),2)=hprima*0.5-(j-1)*(separation_ver2+dv4);

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2)+...
            number_rebars_left(op3),1)=bprima*0.5-(dv*0.5);
    end
end

if number_rebars_inf(op2)>=2
    separation_hor2=round((bprima-((number_rebars_inf(op2))*dv2))/...
        (number_rebars_inf(op2)-1),1);
    
    separation_ver1=round((hprima-((number_rebars_left(op3)+2)*dv3))/...
    (number_rebars_left(op3)+1),1);
    
    % REBAR DISPOSITION________________________________________________
    %__________________________________________________________________
    % rebar coordinates on the inferior part
    for j=1:number_rebars_inf(op2)
        
        disposition_rebar(number_rebars_sup(op1)+j,1)=-bprima*0.5+...
            (dv)*0.5+(j-1)*(separation_hor2+dv2);
        disposition_rebar(number_rebars_sup(op1)+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % left side
    for j=1:number_rebars_left(op3)

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2),2)=...
            hprima*0.5-(dv3)*0.5-(separation_ver1+dv3*0.5)-(j-1)*(separation_ver1+dv3);

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2),1)=...
            -bprima*0.5+(dv*0.5);
    end
else
    separation_ver1=round((hprima-((number_rebars_left(op3)+1)*dv3))/...
    (number_rebars_left(op3)),1);

    % REBAR DISPOSITION________________________________________________
    %__________________________________________________________________
    % rebar coordinates over the inferior part
    for j=1:number_rebars_inf(op2)
        
        disposition_rebar(number_rebars_sup(op1)+j,1)=bprima*0.5-...
            (dv)*0.5;
        disposition_rebar(number_rebars_sup(op1)+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % left side
    for j=1:number_rebars_left(op3)

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2),2)=...
            -hprima*0.5+(dv3)*0.5+(j-1)*(separation_ver1+dv3);

        disposition_rebar(j+number_rebars_sup(op1)+number_rebars_inf(op2),1)=...
            -bprima*0.5+(dv*0.5);

    end

end


