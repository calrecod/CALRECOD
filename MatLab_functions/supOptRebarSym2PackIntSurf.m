function [Mr_col,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
    bestArrangement,bestDisposition,nv4,bestcxy,bestLoad,bestCFA]=...
    supOptRebarSym2PackIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,...
    load_conditions,wac,height,RebarAvailable,ductility,...
    puCostCardBuild,dataCFA,plotRebarDesign)
%-------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
%  bestArrangement,bestDisposition,nv4,bestcxy,bestLoad,bestCFA]=...
%  supOptRebarSym2PackIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,...
%  load_conditions,wac,height,RebarAvailable,ductility,...
%  puCostCardBuild,dataCFA,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design over a 
% rectangular column cross-section in packages of two rebars. Two options 
% available: (1) designs with only one rebar diameter, (2) designs with two
% rebar diameters.
% 
% Note: For the computation of the structural efficiency of each rebar
% design the cross-section is rotated according to each given load combinatio
% so that an interaction diagram is computed for each load.
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
%                               corresponding to the optimal design option
%                               against the mot critical of the given load
%                               conditions
%
%         bestArrangement:      is the list of rebar type of each rebar: 
%                               size [nbars,1] (a number from 1 to 7 by 
%                               default)
%
%         bestDisposition:      is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design option
%
%         bestdiagram:          is the interaction diagram data of the
%                               optimal rebar design 
% 
%         bestcxy:              is the vector containing the neutral axis
%                               depth values for each of the
%                               cross-section's axis, corresponding to the
%                               most critical of the given load conditions
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
%         plotRebarDesign:      is the parameters that indicates if the 
%                               rebar design results are required or not. 
%                               Options are: (1) they are required, 
%                               (2) they are not required
%
%         ductility:            is the parameter that indicates the level
%                               of ductility demand to desing the rebar, 
%                               according to code specifications
%
%         puCostCardBuild:      is a vector containing the parameters
%                               required for the calculation of the unit
%                               cost of a rebar design with a 
%                               "unitCostCardColsRec" 
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
pucb=puCostCardBuild;
pu_col_sym=unitCostCardColsRec(pucb(1),pucb(2),pucb(3),pucb(4),pucb(5),...
                               pucb(6),pucb(7));  
fy=E*0.0021; % yield stress of reinforcing steel 

bp=b-2*rec(1);
hp=h-2*rec(2);

ndiam=length(RebarAvailable(:,1));
noptions=0;
while noptions==0
    bestArea=inf; maxef=1.0;
    for i=1:ndiam % for each type of rebar
        op=i;
        
        ov=RebarAvailable(i,1); % rebar's eight-of-an-inch
        dv=RebarAvailable(i,2); % rebar diameter
        av=(dv)^2*pi/4; % rebar cross-sectional area

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
        % Rebar separation:
        if fdpc<2000 % units (Kg,cm)
            sepMin=max([1.5*2.54, 1.5*dv]); % min rebar separation (cm)
        else         % units (lb,in)
            sepMin=max([1.5, 1.5*dv]); % min rebar separation (in)
        end
        % There is a limit of the number of rebars that can be laid out
        % on each boundary of the cross-section, for each type of rebar
        maxVarillasSup=2*(fix((bp)/(sepMin+2*dv)))+2;
        maxVarillasCos=2*(fix((hp)/(sepMin+2*dv)));
        if 2*maxVarillasSup>nv
            maxVarillasSup=nv/2;
        end
        if 2*maxVarillasCos>(nv-4)
            maxVarillasCos=(nv-4)/2;
        end

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            continue;
        else
            
            bestCost1=nv*av*height*wac*pu_col_sym; % cost of pre 
                                                      % symmetrical design
            
            for type=minVarillasSup:maxVarillasSup
                varSup=type;
                varCos=0.5*(nv-2*varSup);

                % Symmetrical design with only one type of rebar
                arregloVar1=zeros(nv,1)+i;
                [disposicion_varillado]=RebarDisposition2packSym(b,...
                                        h,rec,dv,nv,varCos,varSup);
                if nv~=length(disposicion_varillado)
                    break;
                end
                % Symmetrical design with only one type of rebar
                % with options of two-packs
                [eficiencia,iloadmax,bestLoad1,gamma,bestdiagram1,...
                rotdispositionRebar,rotsection,bestcxy1,cp]=...
                multiDiagAxisColRec(b,h,load_conditions,[i,i,i,i],npdiag,...
                fy,fdpc,beta1,E,varSup,varSup,varCos,varCos,RebarAvailable,...
                disposicion_varillado);
            
                bestEf1=eficiencia(iloadmax,5);
                
                arraySymOriginal=[varSup varSup varCos varCos];
                ab1=varSup*av;
                ab2=varCos*av;
                bestasbar1=ast;
                bestnv1=nv;
                av4_1=[ab1 ab1 ab2 ab2];
                nv4_1=arraySymOriginal;
                bestMr1=eficiencia(iloadmax,4);
                bestDisposition1=disposicion_varillado;
                
                % Symmetrical design with as many as 2 types of rebar
                % with options of two-packs
                
                [av4_2,bestasbar2,bestEf2,bestdiagram2,arregloVar2,...
                bestDisposition2,bestMr2,bestcxy2,bestCost2,bestLoad2,bestCFA2]=...
                sym2typeRebarIntSurf(disposicion_varillado,op,...
                arraySymOriginal,RebarAvailable,b,h,fy,fdpc,beta1,E,...
                load_conditions,wac,height,npdiag,ductility,puCostCardBuild,...
                dataCFA);
            
                bestnv2=nv;
                nv4_2=arraySymOriginal;
                
                % Comparison of best solutions
                if isempty(bestasbar2)==0
                    noptions=noptions+1;
                    if bestasbar2<bestasbar1
                        if bestasbar2<bestArea
                            bestdiagram=bestdiagram2;
                            bestDisposition=bestDisposition2;
                            bestArrangement=arregloVar2;
                            bestLoad=bestLoad2;
                            bestArea=bestasbar2;
                            bestCost=bestCost2;
                            bestCFA=bestCFA2;
                            bestEf=bestEf2;
                            Mr_col=bestMr2;
                            
                            bestnv=bestnv2;
                            nv4=nv4_2;
                            av4=av4_2;
                            bestcxy=bestcxy2;
                        end
                    elseif bestasbar1<=bestasbar2
                        if bestasbar1<bestArea && bestEf1<=maxef
                            bestdiagram=bestdiagram1;
                            bestDisposition=bestDisposition1;
                            bestArrangement=arregloVar1;
                            bestArea=bestasbar1;
                            bestLoad=bestLoad1;
                            bestCost=bestCost1;
                            bestCFA=1;
                            bestEf=bestEf1;
                            Mr_col=bestMr1;
                            
                            bestnv=bestnv1;
                            nv4=nv4_1;
                            av4=av4_1;
                            bestcxy=bestcxy1;
                        end
                    end
                    vx1Ec=nv4(1);
                    vx2Ec=nv4(2); 
                    vy1Ec=nv4(3);
                    vy2Ec=nv4(4);

                    av1Ec=av4(1); 
                    av2Ec=av4(2); 
                    av3Ec=av4(3);
                    av4Ec=av4(4);
                elseif isempty(bestasbar2)==1
                    if bestasbar1<bestArea && bestEf1<=maxef
                        noptions=noptions+1;
                        
                        bestdiagram=bestdiagram1;
                        bestDisposition=bestDisposition1;
                        bestArrangement=arregloVar1;
                        bestArea=bestasbar1;
                        bestLoad=bestLoad1;
                        bestCost=bestCost1;
                        bestCFA=1;
                        bestEf=bestEf1;
                        
                        Mr_col=bestMr1;
                        
                        bestnv=bestnv1;
                        nv4=nv4_1;
                        av4=av4_1;
                        bestcxy=bestcxy1;

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
                
            end
        end
    end
    
    if noptions==0
        fprintf('\nThe dimensions of the columns are too small for rebar.\n');
        fprintf('with this rebar prototype Sym-2pack.\n');
        bestCFA=[];
        bestdiagram=bestdiagram2;
        bestDisposition=[];
        bestArrangement=[];
        bestLoad=[];
        bestArea=[];
        bestCost=[];
        bestEf=[];
        Mr_col=[];
        bestnv=[];
        nv4=[];
        bestcxy=[];
        
        break;
    end
end
% Plotting interaction diagrams and reinforced cross-section (if required)
if plotRebarDesign==1 && noptions>0
    section=[0.5*b 0.5*h;
             -0.5*b 0.5*h;
             -0.5*b -0.5*h;
              0.5*b -0.5*h;
              0.5*b 0.5*h];
          
    PlotRotRecSecRebarCols(section,bestLoad,bestdiagram,...
                    bestDisposition,b,h,bestArrangement);
end