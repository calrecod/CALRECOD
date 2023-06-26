function [Mr_col,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
    bestArrangement,bestDisposition,nv4,bestcxy,bestCP,bestLoad]=...
    Asym4DiamIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,RebarAvailable,...
    height,wac,load_conditions,pu_asym_cols,ductility)
%-------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
%  bestArrangement,bestDisposition,nv4,bestcxy,bestCP,bestLoad]=...
%  Asym4DiamIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,RebarAvailable,...
%  height,wac,load_conditions,pu_asym_cols,ductility)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars asymmetrically 
% over a column cross-section through linear search. The distributions of
% rebars is asymmetrical in quantity of rebars as well as in rebar diameter.
% As many as four rebar diameters are allowed to be simulataneously placed.
% The functions "asym1typeRebarIntSurf" and "asymSym4DiamIntSurf" are used
% to compute such permutations of rebar diameters and numbers of rebars on 
% each of the four cross-section boundaries.
% 
% NOTE: The structural efficiency of each rebar design is determined by
% rotating the cross-section according to each given load combination so
% that their corresponding interaction diagram is computed with respect to
% the action axis of each load condition.
%
% OUTPUT: Mr_col:               are the final resistant bending moment for
%                               both axis directions of the optimal designed 
%                               cross-section
%
%         h:                    modified cross-section height in case it is
%                               modified through the optimization process to
%                               comply the given restrictions of min separation
%                               of rebars
%
%         bestArea:             is the optimal rebar area
%
%         bestCost:             is the cost of the optimal design option
%
%         bestnv:               is the total number of rebars over the 
%                               cross-section corresponding to the optimal
%                               design option
%
%         bestEf:               is the critical structural efficiency 
%                               corresponding to the optimal rebar design
%                               against the most critical of the given load
%                               conditions
%
%         bestArrangement:      is the list of rebar type of each rebar: 
%                               size [nbars,1] (a number from 1 to n#-diam
%
%         bestDisposition:      is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design option
%         bestdiagram:          is the interaction diagram data of the 
%                               optimal rebar design
%
%         bestcxy:              is a vector containing the neutral axis
%                               depth values corresponding to the
%                               most critical load condition for each of
%                               the two cross-section axis
%
%         bestCP:               is a vector containing the Plastic Center
%                               depth values for each of the two
%                               cross-section axis (considering the
%                               asymmetry of the reinforcement)
%
%         bestLoad:             is the resultant most critical load 
%                               combination, in format: [n-load, Pu, Mu],
%                               where Mu = sqrt( Mux^2 + Muy^2 )
%
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
%
%         act:                  optima ISR reinforcement area
%
%         sepMin:               min separation of rebars constraint
%
%         E:                    Elasticity Modulus of reinforcement steel 
%
%         npdiag:               number of points to compute for the 
%                               interaction diagram
%
%         load_conditions:      load conditions for the column cross section:
%                               size = [nload,4] in format [nload,Pu,Mux,Muy]
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to the ACI 318-19 code
%
%         beta1:                is determined as specified by code (see 
%                               Documentation)
%
%         pu_asym_cols:         is the average unit construction cost for
%                               this rebar prototype
%
%         ductility:            is a parameter that indicates which
%                               ductility demand level is required for the
%                               reinforcement designs. A number between 1
%                               to 3 (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
fc=fdpc/0.85;
fy=E*0.0021; % Yield stress of reinforcing steel
pu_col_sym=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93];

bp=b-2*rec(1);
hp=h-2*rec(2);
            
ndiam=length(RebarAvailable(:,1));
noptions=0;
while noptions==0
    bestArea=inf;
    for i=1:ndiam % for each type of rebar
        
        % Completing data-base table:
        op=i;
        
        ov=RebarAvailable(i,1);
        dv=RebarAvailable(i,2);
        av=(dv)^2*pi/4;

        nv=1;
        ast=av*nv;
        while(ast<act)
            nv=nv+1;
            ast=av*nv;
        end
        if (mod(nv,2)~=0)
            nv=nv+1;
            ast=av*nv;
        end
        % Min rebar separation:
        if fc<2000 % units: kg,cm
            sepMin=max([1.5*2.54, 3/2*dv]);
        else       % units: lb,in
            sepMin=max([1.5, 3/2*dv]);
        end
        % There is a limit of the number of rebars that can be laid out
        % on each boundary of the cross-section, for each type of rebar
        maxVarillasSup=fix((bp)/(sepMin+dv))+1;
        maxVarillasCos=fix((hp)/(sepMin+dv))-1;

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            
            continue;
        elseif 2*maxVarillasSup>nv
            continue;
        else
            
            costo=nv*av*height*wac*pu_col_sym(i);
            
            for type=minVarillasSup:maxVarillasSup
                varSup=type;
                varCos=0.5*(nv-2*varSup);

                arregloVar=zeros(nv,1)+i;
                [disposicion_varillado]=RebarDisposition(b,...
                     h,rec,dv,nv,varCos,varSup);
                if nv~=length(disposicion_varillado)
                    break;
                end
                nvxy=[varSup varCos];
                arraySymOriginal=[varSup varSup varCos varCos];
                
                % Asymmetrical design with only one type of rebar
                [av4_1,nv4_1,arregloVar1,bestDisposition1,...
                bestnv1,bestMr1,bestEf1,bestcxy1,bestCP1,bestasbar1,...
                bestdiagram11,bestCost1,maxLoad1]=asym1typeRebarIntSurf...
                (fdpc,nvxy,arraySymOriginal,b,h,rec,RebarAvailable,op,av,...
                npdiag,costo,pu_asym_cols,wac,height,load_conditions,...
                ductility,beta1);
                
                if isempty(nv4_1)==0
                    noptions=noptions+1;
                    % Asymmetrical design with as many as 4 types of rebar
                    % also asymmetrical in number of rebars
                    [av4_2,bestasbar2,bestEf2,bestdiagram21,...
                    arregloVar2,bestDisposition2,bestMr2,...
                    bestcxy2,bestCP2,bestCost2,maxLoad2]=asymSym4DiamIntSurf...
                    (bestDisposition1,op,nv4_1,RebarAvailable,b,h,...
                    fy,fdpc,beta1,E,pu_asym_cols,height,wac,load_conditions,...
                    npdiag,ductility);
                
                    bestnv2=bestnv1;
                    nv4_2=nv4_1;
                else
                    bestasbar2=[];
                end
                
                if isempty(bestasbar2)==0 && isempty(bestasbar1)==0
                    noptions=noptions+1;
                    if bestasbar2<bestasbar1
                        if bestasbar2<bestArea
                            bestdiagram=bestdiagram21;
                            bestLoad=maxLoad2;
                            bestDisposition=bestDisposition2;
                            bestArrangement=arregloVar2;
                            bestArea=bestasbar2;
                            bestCost=bestCost2;
                            bestEf=bestEf2;
                            Mr_col=bestMr2;
                            bestnv=bestnv2;
                            nv4=nv4_2;
                            av4=av4_2;
                            bestcxy=bestcxy2;
                            bestCP=bestCP2;
                            
                            vx1Ec=nv4(1);
                            vx2Ec=nv4(2); 
                            vy1Ec=nv4(3);
                            vy2Ec=nv4(4);

                            av1Ec=av4(1); 
                            av2Ec=av4(2); 
                            av3Ec=av4(3);
                            av4Ec=av4(4);
                        end
                    elseif bestasbar1<=bestasbar2
                        if bestasbar1<bestArea
                            bestdiagram=bestdiagram11;
                            bestLoad=maxLoad1;
                            bestDisposition=bestDisposition1;
                            bestArrangement=arregloVar1;
                            bestArea=bestasbar1;
                            bestCost=bestCost1;
                            bestEf=bestEf1;
                            Mr_col=bestMr1;
                            bestnv=bestnv1;
                            nv4=nv4_1;
                            av4=av4_1;
                            bestcxy=bestcxy1;
                            bestCP=bestCP1;
                            
                            vx1Ec=nv4(1);
                            vx2Ec=nv4(2); 
                            vy1Ec=nv4(3);
                            vy2Ec=nv4(4);

                            av1Ec=av4(1); 
                            av2Ec=av4(2); 
                            av3Ec=av4(3);
                            av4Ec=av4(4);
                        end
                    end
                    
                elseif isempty(bestasbar2)==1 && isempty(bestasbar1)==0
                    noptions=noptions+1;
                    if bestasbar1<bestArea
                        bestdiagram=bestdiagram11;
                        bestLoad=maxLoad1;
                        bestDisposition=bestDisposition1;
                        bestArrangement=arregloVar1;
                        bestArea=bestasbar1;
                        bestCost=bestCost1;
                        bestEf=bestEf1;
                        Mr_col=bestMr1;
                        bestnv=bestnv1;
                        nv4=nv4_1;
                        av4=av4_1;
                        bestcxy=bestcxy1;
                        bestCP=bestCP1;

                        vx1Ec=nv4(1);
                        vx2Ec=nv4(2); 
                        vy1Ec=nv4(3);
                        vy2Ec=nv4(4);

                        av1Ec=av4(1); 
                        av2Ec=av4(2); 
                        av3Ec=av4(3);
                        av4Ec=av4(4);
                    end
                elseif isempty(bestasbar2)==1 && isempty(bestasbar1)==1
                    continue;
                end
            end
        end
    end

    if noptions==0
        fprintf('\nDimensions of column are too small.\n');
        fprintf('The cross-section dimensions should be increased.\n');
        
        bestdiagram=[];
        bestDisposition=[];
        bestArrangement=[];
        bestArea=[];
        bestCost=[];
        bestEf=[];
        Mr_col=[];
        bestnv=[];
        nv4=[];
        bestcxy=[];
        bestCP=[];
        bestLoad=[];
        break;

    end
end