function [dispositionRebar]=RebarDisposition(b,...
    h,rec,dv,nv,varillas_costado,varillas_lecho_alto)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar]=RebarDisposition(b,...
%   h,rec,dv,nv,varillas_costado,varillas_lecho_alto)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local position coordinates of a symmetric rebar
%          design option.
% 
% OUTPUT: dispositionRebar:     are the local position coordinates of the 
%                               symmetric rebar option
%
% INPUT:  b,h:                  given cross-section dimensions
%
%         dv,nv:                are the rebar diameter of the current
%                               option, the number of rebars
%
%         varCos,varSup:        are the number rebars vertically of the 
%                               cross-section (along the cross-section h
%                               height dimension) and the number of rebars
%                               horizontally (along the cross-section b 
%                               width dimension), respectively
%
%         rec:                  is the concrete cover for both cross-section
%                               axis directions: [coverX,coverY]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%------------------------------------------------------------------%%%

bprima=b-2*rec(1);
hprima=h-2*rec(2);

separacion_h=round((hprima-((varillas_costado+2)*dv))/...
    (varillas_costado+1),1);
separacion_b=round((bprima-((varillas_lecho_alto*dv)))/...
    (varillas_lecho_alto-1),1);

dispositionRebar=zeros(nv,2);

% Horizontal rebar disposition_________________________________________
for j=1:varillas_lecho_alto

    % rebar coordinates on the upper boundary
    dispositionRebar(j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
        (separacion_b+dv);
    dispositionRebar(j,2)=hprima*0.5-(dv*0.5);
    
    % rebar coordinates on the lower boundary
    dispositionRebar(varillas_lecho_alto+j,1)=-bprima*0.5+...
        (dv)*0.5+(j-1)*(separacion_b+dv);
    dispositionRebar(varillas_lecho_alto+j,2)=-hprima*0.5+...
        (dv*0.5);

end

% Vertical rebar disposition___________________________________________
for j=1:varillas_costado
    
    % right boundary
    dispositionRebar(j+2*varillas_lecho_alto,2)=hprima*0.5-...
        (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
        (separacion_h+dv);
    dispositionRebar(j+2*varillas_lecho_alto,1)=bprima*0.5-...
        (dv*0.5);
    
    % left boundary
    dispositionRebar(j+2*varillas_lecho_alto+varillas_costado,2)=...
        hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-(j-1)*(separacion_h+dv);
    dispositionRebar(j+2*varillas_lecho_alto+varillas_costado,1)=...
        -bprima*0.5+(dv*0.5);

end
