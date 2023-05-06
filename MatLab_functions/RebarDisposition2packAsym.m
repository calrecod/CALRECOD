function [dispositionRebar]=RebarDisposition2packAsym(b,...
            h,rec,nv,varsup,varsup2,varcos,varcos2,RebarAvailable,...
            op1,op2,op3,op4)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar]=RebarDisposition2packAsym(b,...
%  h,rec,nv,varsup,varsup2,varcos,varcos2,RebarAvailable,...
%  op1,op2,op3,op4)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the local position coordinates of an asymmetrical 
% rebar design over a rectangular reinforced concrete cross-section in 
% packages of two rebars.
% 
% OUTPUT: dispositionRebar:     are the local position coordinates of the 
%                               symmetric rebar option
%
% INPUT:  b,h:                  given cross-section dimensions
%
%         nv:                   the total number of rebars 
%
%         varCos,varCos2:       are the number rebars placed vertically 
%                               over the cross-section (along the 
%                               cross-section h dimension) - left and right
%                               respectively
%
%         varSup,varSup2:       are the number rebars placed horizontally 
%                               over the cross-section (along the 
%                               cross-section b dimension) - upper and 
%                               lower boundary, respectively
%
%         rec:                  is the concrete cover for both cross-section
%                               axis directions: [coverX,coverY]
%
%         RebarAvailable:       Rebar data base. Rebar diameters available,
%                               their respective market numbers and their 
%                               cross-section area, mainly
%
%         op1,op2,op3,op4:      are the rebar diameters' indexes (from the
%                               RebarAvailable database table - rows)
%                               over the cross-section's upper, lower, left
%                               and right boundaries.
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-08-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%------------------------------------------------------------------%%%

bprima=b-2*rec(1);
hprima=h-2*rec(2);

% Rebar separation on the left boundary ------------------------------- 
dv3=RebarAvailable(op3,2);
if varcos>=4 && mod(varcos,2)==0
    separacion_h=round((hprima-((varcos+2)*dv3))/...
        (varcos/2+1),1);
elseif varcos>=4 && mod((varcos-1)/2,2)==0
    separacion_h=round((hprima-((varcos+2)*dv3))/...
        (fix(varcos/2)+2),1);
elseif varcos>=4 && mod((varcos-1)/2,2)==1
    separacion_h=round((hprima-((varcos+2)*dv3))/...
        (fix(varcos/2)+1),1);
else
    separacion_h=round((hprima-((varcos+2)*dv3))/...
        (varcos+1),1);
    
end

% Rebar separation on the right boundary -------------------------------
dv4=RebarAvailable(op4,2);
if varcos2>=4 && mod(varcos2,2)==0
    separacion_h2=round((hprima-((varcos2+2)*dv4))/...
        (varcos2/2+1),1);
elseif varcos2>=4 && mod((varcos2-1)/2,2)==0
    separacion_h2=round((hprima-((varcos2+2)*dv4))/...
        (fix(varcos2/2)+2),1);
elseif varcos2>=4 && mod((varcos2-1)/2,2)==1
    separacion_h2=round((hprima-((varcos2+2)*dv4))/...
        (fix(varcos2/2)+1),1);
else
    separacion_h2=round((hprima-((varcos2+2)*dv4))/...
        (varcos2+1),1);
    
end

% Rebar separation on the upper boundary -------------------------------
dv1=RebarAvailable(op1,2);
if varsup>=4 && mod(varsup,2)==0
    separacion_b=round((bprima-((varsup*dv1)))/...
        (varsup/2-1),1)+dv1;
    
elseif varsup>=4 && mod((varsup-1)/2,2)==0
    separacion_b=round((bprima-((varsup*dv1)))/...
        (fix(varsup/2)),1)+dv1;
elseif varsup>=4 && mod((varsup-1)/2,2)==1
    separacion_b=round((bprima-((varsup*dv1)))/...
        (fix(varsup/2)-1),1)+dv1;
else
    separacion_b=round((bprima-((varsup*dv1)))/...
        (varsup-1),1)+dv1;
end

% Rebar separation on the lower boundary -------------------------------
dv2=RebarAvailable(op2,2);
if varsup2>=4 && mod(varsup2,2)==0
    separacion_b2=round((bprima-((varsup2*dv2)))/...
        (varsup2/2-1),1)+dv2;
    
elseif varsup2>=4 && mod((varsup2-1)/2,2)==0
    separacion_b2=round((bprima-((varsup2*dv2)))/...
        (fix(varsup2/2)),1)+dv2;
elseif varsup2>=4 && mod((varsup2-1)/2,2)==1
    separacion_b2=round((bprima-((varsup2*dv2)))/...
        (fix(varsup2/2)-1),1)+dv2;
else
    separacion_b2=round((bprima-((varsup2*dv2)))/...
        (varsup2-1),1)+dv2;
end

dispositionRebar=zeros(nv,2);

%% Horizontal rebar disposition-----------------------------------------
% Upper boundary -------------------------------------------------------
if varsup>=4
    if mod(varsup,2)==0
        for j=1:2:varsup

            % rebar coordinates on the upper boundary
            dispositionRebar(j,1)=-bprima*0.5+(dv1)*0.5+fix((j-1)/2)*...
                (separacion_b+dv1);
            dispositionRebar(j,2)=hprima*0.5-(dv1*0.5);

            dispositionRebar(j+1,1)=-bprima*0.5+(dv1)*0.5+dv1+fix((j-1)/2)*...
                (separacion_b+dv1);
            dispositionRebar(j+1,2)=hprima*0.5-(dv1*0.5);
        end
    elseif mod(varsup,2)==1 && ... % two-packs and one single 
           mod((varsup-1)/2,2)==0  % rebar at the middle
                                     
        npacks=fix(varsup/2)+1;
        for j=1:npacks

            if j~=fix(npacks/2)+1 % for the rebars in two-packs at the left 
                                  % and right side
                if j<fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv1*0.5);

                    dispositionRebar(2*j,1)=-bprima*0.5+(dv1)*0.5+dv1+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv1*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-2,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j-2,2)=hprima*0.5-(dv1*0.5);

                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv1)*0.5+dv1+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv1*0.5);
                end
            else % for the single rebar at the middle

                % rebar coordinates on the upper boundary
                dispositionRebar(2*j-1,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                    (separacion_b+dv1);
                dispositionRebar(2*j-1,2)=hprima*0.5-(dv1*0.5);
            end

        end
    elseif mod(varsup,2)==1 && mod((varsup-1)/2,2)==1  % two-packs and one
                                                       % and on three-pack
        npacks=fix(varsup/2);
        for j=1:npacks
            if j~=fix(npacks/2)+1
                if j<fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv1*0.5);

                    dispositionRebar(2*j,1)=-bprima*0.5+(dv1)*0.5+dv1+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv1*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv1*0.5);

                    dispositionRebar(2*j+1,1)=-bprima*0.5+(dv1)*0.5+dv1+(j-1)*...
                        (separacion_b+dv1);
                    dispositionRebar(2*j+1,2)=hprima*0.5-(dv1*0.5);
                end
            else % for the three-pack rebar at the middle

                % rebar coordinates on the upper boundary
                dispositionRebar(2*j-1,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
                    (separacion_b+dv1);
                dispositionRebar(2*j-1,2)=hprima*0.5-(dv1*0.5);

                dispositionRebar(2*j,1)=-bprima*0.5+(dv1)*0.5+dv1+(j-1)*...
                    (separacion_b+dv1);
                dispositionRebar(2*j,2)=hprima*0.5-(dv1*0.5);
                
                dispositionRebar(2*j+1,1)=-bprima*0.5+(dv1)*0.5+2*dv1+(j-1)*...
                    (separacion_b+dv1);
                dispositionRebar(2*j+1,2)=hprima*0.5-(dv1*0.5);
            end
        end
    end
else
    for j=1:varsup
        % rebar coordinates on the upper boundary
        dispositionRebar(j,1)=-bprima*0.5+(dv1)*0.5+(j-1)*...
            (separacion_b+dv1);
        dispositionRebar(j,2)=hprima*0.5-(dv1*0.5);
    end
end


% Lower boundary -------------------------------------------------------
if varsup2>=4
    if mod(varsup2,2)==0
        for j=1:2:varsup2

            % rebar coordinates on the lower boundary
            dispositionRebar(varsup+j,1)=-bprima*0.5+...
                (dv2)*0.5+fix((j-1)/2)*(separacion_b2+dv2);
            dispositionRebar(varsup+j,2)=-hprima*0.5+...
                (dv2*0.5);

            dispositionRebar(varsup+j+1,1)=-bprima*0.5+...
                (dv2)*0.5+dv2+fix((j-1)/2)*(separacion_b2+dv2);
            dispositionRebar(varsup+j+1,2)=-hprima*0.5+...
                (dv2*0.5);
        end
    elseif mod(varsup2,2)==1 && ... % two-packs and one single 
           mod((varsup2-1)/2,2)==0  % rebar at the middle
                                     
        npacks=fix(varsup2/2)+1;
        for j=1:npacks
            if j~=fix(npacks/2)+1 % for the rebars in two-packs at the left 
                                  % and right side
                if j<fix(npacks/2)+1
                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                    (dv2)*0.5+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv2*0.5);

                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                        (dv2)*0.5+dv2+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv2*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-2,1)=-bprima*0.5+...
                    (dv2)*0.5+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j-2,2)=-hprima*0.5+...
                        (dv2*0.5);

                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                        (dv2)*0.5+dv2+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv2*0.5);
                end
            else % for the single rebar at the middle

                % rebar coordinates on the lower boundary
                dispositionRebar(2*j-1+varsup,1)=-bprima*0.5+(dv2)*0.5+(j-1)*...
                    (separacion_b2+dv2);
                dispositionRebar(2*j-1+varsup,2)=-hprima*0.5+(dv2*0.5);
            end

        end
    elseif mod(varsup2,2)==1 && mod((varsup2-1)/2,2)==1  % two-packs and one
                                                       % and on three-pack
        npacks=fix(varsup2/2);
        for j=1:npacks
            if j~=fix(npacks/2)+1
                if j<fix(npacks/2)+1
                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                    (dv2)*0.5+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv2*0.5);

                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                        (dv2)*0.5+dv2+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv2*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                    (dv2)*0.5+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv2*0.5);

                    dispositionRebar(varsup+2*j+1,1)=-bprima*0.5+...
                        (dv2)*0.5+dv2+(j-1)*(separacion_b2+dv2);
                    dispositionRebar(varsup+2*j+1,2)=-hprima*0.5+...
                        (dv2*0.5);
                end
            else % for the three-pack rebar at the middle

                % rebar coordinates on the lower boundary
                dispositionRebar(2*j-1+varsup,1)=-bprima*0.5+(dv2)*0.5+(j-1)*...
                    (separacion_b2+dv2);
                dispositionRebar(2*j-1+varsup,2)=-hprima*0.5+(dv2*0.5);
                
                dispositionRebar(2*j+varsup,1)=-bprima*0.5+(dv2)*0.5+dv2+(j-1)*...
                    (separacion_b2+dv2);
                dispositionRebar(2*j+varsup,2)=-hprima*0.5+(dv2*0.5);
                
                dispositionRebar(2*j+1+varsup,1)=-bprima*0.5+(dv2)*0.5+2*...
                    dv2+(j-1)*(separacion_b2+dv2);
                dispositionRebar(2*j+1+varsup,2)=-hprima*0.5+(dv2*0.5);
            
            end
        end
    end
else
    
    for j=1:varsup2
        % rebar coordinates on the lower boundary
        dispositionRebar(j+varsup,1)=-bprima*0.5+(dv2)*0.5+(j-1)*...
            (separacion_b2+dv2);
        dispositionRebar(j+varsup,2)=-hprima*0.5+(dv2*0.5);
    end
end

%% Vertical rebar disposition--------------------------------------------
% Left boundary---------------------------------------------------------
if mod(varcos,2)==0
    for j=1:2:varcos

        % left boundary
        dispositionRebar(j+varsup+varsup2,2)=...
            hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-fix((j-1)/2)*...
            (separacion_h+dv3);
        dispositionRebar(j+varsup+varsup2,1)=...
            -bprima*0.5+(dv3*0.5);

        dispositionRebar(j+1+varsup+varsup2,2)=...
            hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-fix((j-1)/2)*...
            (separacion_h+dv3);
        dispositionRebar(j+1+varsup+varsup2,1)=...
            -bprima*0.5+(dv3*0.5);

    end
else % two-packs and one single rebar at the
     % middle
    if  varcos>4
        if mod((varcos-1)/2,2)==0
            npacks=fix(varcos/2)+1;
            for j=1:npacks
                % for the rebars in two-packs at the upper and lower side
                if j~=fix(npacks/2)+1
                    if j<fix(npacks/2)+1
                        % left boundary
                        dispositionRebar(2*j-1+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j-1+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);

                        dispositionRebar(2*j+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);
                    elseif j>fix(npacks/2)+1
                        % left boundary
                        dispositionRebar(2*j-2+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j-2+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);

                        dispositionRebar(2*j-1+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j-1+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);
                    end
                else
                    % left boundary
                    dispositionRebar(2*j-1+varsup+varsup2,2)=...
                        hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-(j-1)*...
                        (separacion_h+dv3);
                    dispositionRebar(2*j-1+varsup+varsup2,1)=...
                        -bprima*0.5+(dv3*0.5);
                end
            end
        elseif mod((varcos-1)/2,2)==1
            npacks=fix(varcos/2);
            for j=1:npacks
                % for the rebars in two-packs at the upper and lower side
                if j~=fix(npacks/2)+1
                    if j<fix(npacks/2)+1
                        % left boundary
                        dispositionRebar(2*j-1+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j-1+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);

                        dispositionRebar(2*j+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);
                    elseif j>fix(npacks/2)+1
                        % left boundary
                        dispositionRebar(2*j+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);

                        dispositionRebar(2*j+1+varsup+varsup2,2)=...
                            hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-...
                            (j-1)*(separacion_h+dv3);
                        dispositionRebar(2*j+1+varsup+varsup2,1)=...
                            -bprima*0.5+(dv3*0.5);
                    end
                else
                    % left boundary
                    dispositionRebar(2*j-1+varsup+varsup2,2)=...
                        hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-(j-1)*...
                        (separacion_h+dv3);
                    dispositionRebar(2*j-1+varsup+varsup2,1)=...
                        -bprima*0.5+(dv3*0.5);
                    
                    dispositionRebar(2*j+varsup+varsup2,2)=...
                        hprima*0.5-(dv3)*0.5-dv3-(separacion_h+dv3*0.5)-...
                        (j-1)*(separacion_h+dv3);
                    dispositionRebar(2*j+varsup+varsup2,1)=...
                        -bprima*0.5+(dv3*0.5);
                    
                    dispositionRebar(2*j+1+varsup+varsup2,2)=...
                        hprima*0.5-(dv3)*0.5-2*dv3-(separacion_h+dv3*0.5)-...
                        (j-1)*(separacion_h+dv3);
                    dispositionRebar(2*j+1+varsup+varsup2,1)=...
                        -bprima*0.5+(dv3*0.5);
                end
            end
        end
    else % if the number of rebars in vertical are less than four
        for j=1:varcos
            % left boundary
            dispositionRebar(j+varsup+varsup2,2)=...
                hprima*0.5-(dv3)*0.5-(separacion_h+dv3*0.5)-(j-1)*...
                (separacion_h+dv3);
            dispositionRebar(j+varsup+varsup2,1)=...
                -bprima*0.5+(dv3*0.5);
        end
    end
end

% Right boundary--------------------------------------------------------
if mod(varcos2,2)==0
    for j=1:2:varcos2
        
        % right boundary
        dispositionRebar(j+varsup+varsup2+varcos,2)=hprima*0.5-...
            (dv4)*0.5-(separacion_h2+dv4*0.5)-fix((j-1)/2)*...
            (separacion_h2+dv4);
        dispositionRebar(j+varsup+varsup2+varcos,1)=bprima*0.5-...
            (dv4*0.5);

        dispositionRebar(j+1+varsup+varsup2+varcos,2)=hprima*0.5-...
            (dv4)*0.5-dv4-(separacion_h2+dv4*0.5)-fix((j-1)/2)*...
            (separacion_h2+dv4);
        dispositionRebar(j+1+varsup+varsup2+varcos,1)=bprima*0.5-...
            (dv4*0.5);
    end
else % two-packs and one single rebar at the
     % middle
    if  varcos2>4
        if mod((varcos2-1)/2,2)==0
            npacks=fix(varcos2/2)+1;
            for j=1:npacks
                % for the rebars in two-packs at the upper and lower side
                if j~=fix(npacks/2)+1
                    if j<fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j-1+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-...
                            (j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j-1+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);

                        dispositionRebar(2*j+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-dv4-(separacion_h2+...
                            dv4*0.5)-(j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);
                    elseif j>fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j-2+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-...
                            (j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j-2+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);

                        dispositionRebar(2*j-1+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-dv4-(separacion_h2+dv4*...
                            0.5)-(j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j-1+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);
                    end
                else
                    % right boundary
                    dispositionRebar(2*j-1+varsup+varsup2+varcos,2)=...
                        hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-(j-1)*...
                        (separacion_h2+dv4);
                    dispositionRebar(2*j-1+varsup+varsup2+varcos,1)=...
                        bprima*0.5-(dv4*0.5);
                end
            end
        elseif mod((varcos2-1)/2,2)==1
            npacks=fix(varcos2/2);
            for j=1:npacks
                % for the rebars in two-packs at the upper and lower side
                if j~=fix(npacks/2)+1
                    if j<fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j-1+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-...
                            (j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j-1+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);

                        dispositionRebar(2*j+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-dv4-(separacion_h2+dv4*...
                            0.5)-(j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);
                    elseif j>fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-...
                            (j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);

                        dispositionRebar(2*j+1+varsup+varsup2+varcos,2)=...
                            hprima*0.5-(dv4)*0.5-dv4-(separacion_h2+dv4*...
                            0.5)-(j-1)*(separacion_h2+dv4);
                        dispositionRebar(2*j+1+varsup+varsup2+varcos,1)=...
                            bprima*0.5-(dv4*0.5);
                    end
                else
                    % right boundary
                    dispositionRebar(2*j-1+varsup+varsup2+varcos,2)=...
                        hprima*0.5-(dv4)*0.5-(separacion_h2+dv4*0.5)-(j-1)*...
                        (separacion_h2+dv4);
                    dispositionRebar(2*j-1+varsup+varsup2+varcos,1)=...
                        bprima*0.5-(dv4*0.5);
                    
                    dispositionRebar(2*j+varsup+varsup2+varcos,2)=hprima*...
                        0.5-(dv4)*0.5-dv4-(separacion_h2+dv4*0.5)-(j-1)*...
                        (separacion_h2+dv4);
                    dispositionRebar(2*j+varsup+varsup2+varcos,1)=bprima*...
                        0.5-(dv4*0.5);
                    
                    dispositionRebar(2*j+1+varsup+varsup2+varcos,2)=hprima*...
                        0.5-(dv4)*0.5-2*dv4-(separacion_h2+dv4*0.5)-(j-1)*...
                        (separacion_h2+dv4);
                    dispositionRebar(2*j+1+varsup+varsup2+varcos,1)=bprima*...
                        0.5-(dv4*0.5);
                end
            end
        end
    else % if the number of rebars in vertical are less than four
        for j=1:varcos2
            % right boundary
            dispositionRebar(j+varsup+varsup2+varcos,2)=hprima*0.5-...
                (dv4)*0.5-(separacion_h2+dv4*0.5)-(j-1)*...
                (separacion_h2+dv4);
            dispositionRebar(j+varsup+varsup2+varcos,1)=bprima*0.5-...
                (dv4*0.5);
        end
    end
end