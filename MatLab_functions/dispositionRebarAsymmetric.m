function [dispositionRebar,sephor1,sephor2,...
    sepver1,sepver2]=dispositionRebarAsymmetric(b,...
    h,rec,nv,nRebarsSup,nRebarsInf,nRebarsLeft,...
    nRebarsRight,RebarAvailable,op1,op2,op3,op4)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar,sephor1,sephor2,...
%  sepver1,sepver2]=dispositionRebarAsymmetric(b,...
%  h,rec,nv,nRebarsSup,nRebarsInf,nRebarsLeft,...
%  nRebarsRight,RebarAvailable,op1,op2,op3,op4)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local coordinates of an asymmetric rebar option.
% 
% OUTPUT: dispositionRebar:     are the local coordinates of the optimal 
%                               rebar option
%
%         sephor1,
%         sephor2,
%         sepver1,
%         sepver2:              resultant rebar separation to be compared 
%                               with the minimum one (upper, lower, left
%                               right boundary), respectively
%
% INPUT:  b,h:                  cross-section initial dimensions
%
%         rec:                  are the concrete cover values for each axis
%                               direction of the cross-section
%
%         RebarAvailable:       rebar database consisting of an array of 
%                               size: n# x 2, by default in format: 
%                               [#rebar, diam]
%
%         nRebarsSup,
%         nRebarsInf,
%         nRebarsLeft,
%         nRebarsRight:         number of rebars to placed on each of the 
%                               cross-section boundaries
%
%         op1,op2,op3,op4:      resultant types of rebar for each of the 
%                               four cross-section boundaries (upper 
%                               boundary, lower boundary, left side and 
%                               right side, respectively)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

bprima=b-2*rec(1);
hprima=h-2*rec(2);

dv1=RebarAvailable(op1,2);
dv2=RebarAvailable(op2,2);
dv3=RebarAvailable(op3,2);
dv4=RebarAvailable(op4,2);

ndiam=length(RebarAvailable(:,1));
dispositionRebar=zeros(nv,2);
dv=RebarAvailable(ndiam,2);
if nRebarsSup>=2
    
    sephor1=round((bprima-((nRebarsSup)*dv1))/...
        (nRebarsSup-1),1);
    
    sepver2=round(((hprima)-((nRebarsRight+2)*dv4))/...
    (nRebarsRight+1),1);

    % REBAR DISTRIBUTION-----------------------------------------------
    % rebars disposition horizontally (superior)
    for j=1:nRebarsSup
        dispositionRebar(j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
            (sephor1+dv1);
        dispositionRebar(j,2)=hprima*0.5-(dv*0.25);
    end
    
    % rebar disposition in the vertical right part
    for j=1:nRebarsRight
        dispositionRebar(j+nRebarsSup+nRebarsInf+...
            nRebarsLeft,2)=hprima*0.5-(dv4)*0.5-(sepver2+dv4*0.5)-(j-1)*(sepver2+dv4);

        dispositionRebar(j+nRebarsSup+nRebarsInf+...
            nRebarsLeft,1)=bprima*0.5-(dv*0.5);
    end

else
    sepver2=round((hprima-((nRebarsRight+1)*dv4))/...
    (nRebarsRight),1);

    % REBAR DISTRIBUTION-----------------------------------------------
    % rebar disposition horizontally (superior) only one rebar
    for j=1:nRebarsSup
        dispositionRebar(j,1)=-bprima*0.5+(dv)*0.5;
        dispositionRebar(j,2)=hprima*0.5-(dv*0.5);
    end
    
    % rebar disposition over right side
    for j=1:nRebarsRight
        dispositionRebar(j+nRebarsSup+nRebarsInf+...
            nRebarsLeft,2)=hprima*0.5-(j-1)*(sepver2+dv4);

        dispositionRebar(j+nRebarsSup+nRebarsInf+...
            nRebarsLeft,1)=bprima*0.5-(dv*0.5);
    end
end

if nRebarsInf>=2
    sephor2=round((bprima-((nRebarsInf)*dv2))/...
        (nRebarsInf-1),1);
    
    sepver1=round((hprima-((nRebarsLeft+2)*dv3))/...
    (nRebarsLeft+1),1);
    
    % REBAR DISTRIBUTION-----------------------------------------------
    % rebar coordinates on the inferior part
    for j=1:nRebarsInf
        
        dispositionRebar(nRebarsSup+j,1)=-bprima*0.5+...
            (dv)*0.5+(j-1)*(sephor2+dv2);
        dispositionRebar(nRebarsSup+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % rebar coordinates on the left side
    for j=1:nRebarsLeft

        dispositionRebar(j+nRebarsSup+nRebarsInf,2)=...
            hprima*0.5-(dv3)*0.5-(sepver1+dv3*0.5)-(j-1)*(sepver1+dv3);

        dispositionRebar(j+nRebarsSup+nRebarsInf,1)=...
            -bprima*0.5+(dv*0.5);
    end
else
    sepver1=round((hprima-((nRebarsLeft+1)*dv3))/...
    (nRebarsLeft),1);

    % REBAR DISTRIBUTION-----------------------------------------------
    % rebar coordinates over the inferior part
    for j=1:nRebarsInf
        
        dispositionRebar(nRebarsSup+j,1)=bprima*0.5-...
            (dv)*0.5;
        dispositionRebar(nRebarsSup+j,2)=-hprima*0.5+...
            (dv*0.25);

    end
    
    % rebar coordinates over the left side
    for j=1:nRebarsLeft

        dispositionRebar(j+nRebarsSup+nRebarsInf,2)=...
            -hprima*0.5+(dv3)*0.5+(j-1)*(sepver1+dv3);

        dispositionRebar(j+nRebarsSup+nRebarsInf,1)=...
            -bprima*0.5+(dv*0.5);

    end

end