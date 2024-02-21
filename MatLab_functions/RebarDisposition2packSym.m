function [dispositionRebar]=RebarDisposition2packSym(b,...
            h,rec,dv,nv,varcos,varsup)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar]=RebarDisposition2packSym(b,...
%   h,rec,dv,nv,varcos,varsup)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
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
%         varCos,varSup:        are the number of rebars vertically of the 
%                               cross-section (along the cross-section h
%                               height dimension) and the number of rebars
%                               horizontally (along the cross-section b 
%                               width dimension), respectively
%
%         rec:                  is the concrete cover for both cross-section
%                               axis directions: [coverX,coverY]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

bprima=b-2*rec(1);
hprima=h-2*rec(2);

if varcos>=4 && mod(varcos,2)==0
    separacion_h=round((hprima-((varcos+2)*dv))/...
        (varcos/2+1),1);
elseif varcos>=4 && mod((varcos-1)/2,2)==0
    separacion_h=round((hprima-((varcos+2)*dv))/...
        (fix(varcos/2)+2),1);
elseif varcos>=4 && mod((varcos-1)/2,2)==1
    separacion_h=round((hprima-((varcos+2)*dv))/...
        (fix(varcos/2)+1),1);
else
    separacion_h=round((hprima-((varcos+2)*dv))/...
        (varcos+1),1);
end
if varsup>=4 && mod(varsup,2)==0
    separacion_b=round((bprima-((varsup*dv)))/...
        (varsup/2-1),1)+dv;
    
elseif varsup>=4 && mod((varsup-1)/2,2)==0
    separacion_b=round((bprima-((varsup*dv)))/...
        (fix(varsup/2)),1)+dv;
elseif varsup>=4 && mod((varsup-1)/2,2)==1
    separacion_b=round((bprima-((varsup*dv)))/...
        (fix(varsup/2)-1),1)+dv;
else
    separacion_b=round((bprima-((varsup*dv)))/...
        (varsup-1),1);
end
dispositionRebar=zeros(nv,2);

%% Horizontal rebar disposition------------------------------------------
if varsup>=4
    if mod(varsup,2)==0
        for j=1:2:varsup
            % rebar coordinates on the upper boundary
            dispositionRebar(j,1)=-bprima*0.5+(dv)*0.5+fix((j-1)/2)*...
                (separacion_b+dv);
            dispositionRebar(j,2)=hprima*0.5-(dv*0.5);

            dispositionRebar(j+1,1)=-bprima*0.5+(dv)*0.5+dv+fix((j-1)/2)*...
                (separacion_b+dv);
            dispositionRebar(j+1,2)=hprima*0.5-(dv*0.5);

            % rebar coordinates on the lower boundary
            dispositionRebar(varsup+j,1)=-bprima*0.5+...
                (dv)*0.5+fix((j-1)/2)*(separacion_b+dv);
            dispositionRebar(varsup+j,2)=-hprima*0.5+...
                (dv*0.5);

            dispositionRebar(varsup+j+1,1)=-bprima*0.5+...
                (dv)*0.5+dv+fix((j-1)/2)*(separacion_b+dv);
            dispositionRebar(varsup+j+1,2)=-hprima*0.5+...
                (dv*0.5);
        end
    elseif mod(varsup,2)==1 && ... % two-packs and one single 
           mod((varsup-1)/2,2)==0  % rebar at the middle
                      
        npacks=fix(varsup/2)+1;
        for j=1:npacks

            if j~=fix(npacks/2)+1 % for the rebars in two-packs at the left 
                                  % and right side
                if j<fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv*0.5);

                    dispositionRebar(2*j,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv*0.5);

                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                    (dv)*0.5+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv*0.5);

                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                        (dv)*0.5+dv+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-2,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j-2,2)=hprima*0.5-(dv*0.5);

                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv*0.5);

                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-2,1)=-bprima*0.5+...
                    (dv)*0.5+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j-2,2)=-hprima*0.5+...
                        (dv*0.5);

                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                        (dv)*0.5+dv+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv*0.5);
                end
            else % for the single rebar at the middle

                % rebar coordinates on the upper boundary
                dispositionRebar(2*j-1,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j-1,2)=hprima*0.5-(dv*0.5);

                % rebar coordinates on the lower boundary
                dispositionRebar(2*j-1+varsup,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j-1+varsup,2)=-hprima*0.5+(dv*0.5);
            end

        end
    elseif mod(varsup,2)==1 && mod((varsup-1)/2,2)==1  % two-packs and one
                                                       % on three-pack

        npacks=fix(varsup/2);
        for j=1:npacks
            if j~=fix(npacks/2)+1
                if j<fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j-1,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j-1,2)=hprima*0.5-(dv*0.5);

                    dispositionRebar(2*j,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv*0.5);

                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j-1,1)=-bprima*0.5+...
                    (dv)*0.5+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j-1,2)=-hprima*0.5+...
                        (dv*0.5);

                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                        (dv)*0.5+dv+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv*0.5);
                elseif j>fix(npacks/2)+1
                    % rebar coordinates on the upper boundary
                    dispositionRebar(2*j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j,2)=hprima*0.5-(dv*0.5);

                    dispositionRebar(2*j+1,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                        (separacion_b+dv);
                    dispositionRebar(2*j+1,2)=hprima*0.5-(dv*0.5);

                    % rebar coordinates on the lower boundary
                    dispositionRebar(varsup+2*j,1)=-bprima*0.5+...
                    (dv)*0.5+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j,2)=-hprima*0.5+...
                        (dv*0.5);

                    dispositionRebar(varsup+2*j+1,1)=-bprima*0.5+...
                        (dv)*0.5+dv+(j-1)*(separacion_b+dv);
                    dispositionRebar(varsup+2*j+1,2)=-hprima*0.5+...
                        (dv*0.5);
                end
            else % for the three-pack rebar at the middle

                % rebar coordinates on the upper boundary
                dispositionRebar(2*j-1,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j-1,2)=hprima*0.5-(dv*0.5);

                dispositionRebar(2*j,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j,2)=hprima*0.5-(dv*0.5);
                
                dispositionRebar(2*j+1,1)=-bprima*0.5+(dv)*0.5+2*dv+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j+1,2)=hprima*0.5-(dv*0.5);
                
                % rebar coordinates on the lower boundary
                dispositionRebar(2*j-1+varsup,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j-1+varsup,2)=-hprima*0.5+(dv*0.5);
                
                dispositionRebar(2*j+varsup,1)=-bprima*0.5+(dv)*0.5+dv+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j+varsup,2)=-hprima*0.5+(dv*0.5);
                
                dispositionRebar(2*j+1+varsup,1)=-bprima*0.5+(dv)*0.5+2*dv+(j-1)*...
                    (separacion_b+dv);
                dispositionRebar(2*j+1+varsup,2)=-hprima*0.5+(dv*0.5);
            
            end
        end
    end
else
    for j=1:varsup
        % rebar coordinates on the upper boundary
        dispositionRebar(j,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
            (separacion_b+dv);
        dispositionRebar(j,2)=hprima*0.5-(dv*0.5);

        % rebar coordinates on the lower boundary
        dispositionRebar(j+varsup,1)=-bprima*0.5+(dv)*0.5+(j-1)*...
            (separacion_b+dv);
        dispositionRebar(j+varsup,2)=-hprima*0.5+(dv*0.5);
    end
end
%% Vertical rebar disposition--------------------------------------------
if mod(varcos,2)==0
    for j=1:2:varcos

        % right boundary
        dispositionRebar(j+2*varsup,2)=hprima*0.5-...
            (dv)*0.5-(separacion_h+dv*0.5)-fix((j-1)/2)*...
            (separacion_h+dv);
        dispositionRebar(j+2*varsup,1)=bprima*0.5-...
            (dv*0.5);

        dispositionRebar(j+1+2*varsup,2)=hprima*0.5-...
            (dv)*0.5-dv-(separacion_h+dv*0.5)-fix((j-1)/2)*...
            (separacion_h+dv);
        dispositionRebar(j+1+2*varsup,1)=bprima*0.5-...
            (dv*0.5);

        % left boundary
        dispositionRebar(j+2*varsup+varcos,2)=...
            hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-fix((j-1)/2)*...
            (separacion_h+dv);
        dispositionRebar(j+2*varsup+varcos,1)=...
            -bprima*0.5+(dv*0.5);

        dispositionRebar(j+1+2*varsup+varcos,2)=...
            hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-fix((j-1)/2)*...
            (separacion_h+dv);
        dispositionRebar(j+1+2*varsup+varcos,1)=...
            -bprima*0.5+(dv*0.5);

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
                        % right boundary
                        dispositionRebar(2*j-1+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        dispositionRebar(2*j+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-dv-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        % left boundary
                        dispositionRebar(2*j-1+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);

                        dispositionRebar(2*j+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);
                    elseif j>fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j-2+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j-2+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        dispositionRebar(2*j-1+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-dv-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        % left boundary
                        dispositionRebar(2*j-2+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j-2+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);

                        dispositionRebar(2*j-1+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);
                    end
                else
                    % right boundary
                    dispositionRebar(2*j-1+2*varsup,2)=hprima*0.5-...
                        (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j-1+2*varsup,1)=bprima*0.5-...
                        (dv*0.5);

                    % left boundary
                    dispositionRebar(2*j-1+2*varsup+varcos,2)=...
                        hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j-1+2*varsup+varcos,1)=...
                        -bprima*0.5+(dv*0.5);
                end
            end
        elseif mod((varcos-1)/2,2)==1
            npacks=fix(varcos/2);
            for j=1:npacks
                % for the rebars in two-packs at the upper and lower side
                if j~=fix(npacks/2)+1
                    if j<fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j-1+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        dispositionRebar(2*j+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-dv-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        % left boundary
                        dispositionRebar(2*j-1+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j-1+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);

                        dispositionRebar(2*j+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);
                    elseif j>fix(npacks/2)+1
                        % right boundary
                        dispositionRebar(2*j+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        dispositionRebar(2*j+1+2*varsup,2)=hprima*0.5-...
                            (dv)*0.5-dv-(separacion_h+dv*0.5)-(j-1)*...
                            (separacion_h+dv);
                        dispositionRebar(2*j+1+2*varsup,1)=bprima*0.5-...
                            (dv*0.5);

                        % left boundary
                        dispositionRebar(2*j+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);

                        dispositionRebar(2*j+1+2*varsup+varcos,2)=...
                            hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-...
                            (j-1)*(separacion_h+dv);
                        dispositionRebar(2*j+1+2*varsup+varcos,1)=...
                            -bprima*0.5+(dv*0.5);
                    end
                else
                    % right boundary
                    dispositionRebar(2*j-1+2*varsup,2)=hprima*0.5-...
                        (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j-1+2*varsup,1)=bprima*0.5-...
                        (dv*0.5);
                    
                    dispositionRebar(2*j+2*varsup,2)=hprima*0.5-...
                        (dv)*0.5-dv-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j+2*varsup,1)=bprima*0.5-...
                        (dv*0.5);
                    
                    dispositionRebar(2*j+1+2*varsup,2)=hprima*0.5-...
                        (dv)*0.5-2*dv-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j+1+2*varsup,1)=bprima*0.5-...
                        (dv*0.5);

                    % left boundary
                    dispositionRebar(2*j-1+2*varsup+varcos,2)=...
                        hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                        (separacion_h+dv);
                    dispositionRebar(2*j-1+2*varsup+varcos,1)=...
                        -bprima*0.5+(dv*0.5);
                    
                    dispositionRebar(2*j+2*varsup+varcos,2)=...
                        hprima*0.5-(dv)*0.5-dv-(separacion_h+dv*0.5)-...
                        (j-1)*(separacion_h+dv);
                    dispositionRebar(2*j+2*varsup+varcos,1)=...
                        -bprima*0.5+(dv*0.5);
                    
                    dispositionRebar(2*j+1+2*varsup+varcos,2)=...
                        hprima*0.5-(dv)*0.5-2*dv-(separacion_h+dv*0.5)-...
                        (j-1)*(separacion_h+dv);
                    dispositionRebar(2*j+1+2*varsup+varcos,1)=...
                        -bprima*0.5+(dv*0.5);
                end
            end
        end
    else % if the number of rebars in vertical are less than four
        for j=1:varcos
            % right boundary
            dispositionRebar(j+2*varsup,2)=hprima*0.5-...
                (dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                (separacion_h+dv);
            dispositionRebar(j+2*varsup,1)=bprima*0.5-...
                (dv*0.5);

            % left boundary
            dispositionRebar(j+2*varsup+varcos,2)=...
                hprima*0.5-(dv)*0.5-(separacion_h+dv*0.5)-(j-1)*...
                (separacion_h+dv);
            dispositionRebar(j+2*varsup+varcos,1)=...
                -bprima*0.5+(dv*0.5);
        end
    end
end