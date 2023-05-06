function [Mr_col,h,Inertia_xy_modif,bestArea,bestCost,bestdiagram,bestnv,...
    bestEf,bestArrangement,bestDisposition,nv4,bestcxy]=superOptimalRebarSym...
    (b,h,rec,act,E,npdiag,fdpc,beta1,pu_col_sym,load_conditions,...
    wac,height,condition_cracking,ductility,RebarAvailable,plotRebarDesign)

%-------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,Inertia_xy_modif,bestArea,bestCost,bestdiagram,bestnv,...
%  bestEf,bestArrangement,bestDisposition,nv4,bestcxy]=superOptimalRebarSym...
%  (b,h,rec,act,E,npdiag,fdpc,beta1,pu_col_sym,load_conditions,...
%  wac,height,condition_cracking,ductility,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal rebar design symmetrically distributed
% over a rectangular column cross-section. Two options are available:
% (1) only one rebar diameter, (2) two rebar diameters.
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
%         Inertia_xy_modif:     momentum of inertia of the bar reinforced
%                               cross-section for both axis directions, by 
%                               the computation of the cracking mechanisms 
%                               according to the parameter condition_cracking
%
%         bestArea:             is the optimal rebar area
%
%         bestCost:             is the cost of the optimal design option
%
%         bestEf:               is the critical structural efficiency 
%                               corresponding to the optimal design against
%                               the most critical of the given load
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
%         bestnv:               is the total number of rebars corresponding
%                               to the optimal design
%
%         bestdiagram:          is the interaction diagram data of the
%                               optimal rebar design 
% 
%         bestcxy:              is the vector containing the neutral axis
%                               depth values for each of the
%                               cross-section's axis, corresponding to the
%                               most critical of the given load conditions
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
%         pu_col:               is the database of reinforcement assembly
%                               and construction unit cost: format by
%                               default:
%    -----------------------------------------------------------------
%    pu_col=[PU{#4}, PU{#5}, PU{#6}, PU{#8}, PU{#9}, PU{#10}, PU{#12}]
%    -----------------------------------------------------------------
%
%         condition_cracking:   parameter that indicates which cross-section
%                               cracking mechanism will be consider, either 
%                               Cracked or Non-cracked. If the condition 
%                               Non-cracked is set, then the cracking 
%                               mechanism will be neglected by all means
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
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-21
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fy=E*0.0021; % yield stress of reinforcing steel
puSym2cols=1.1*sum(pu_col_sym)/length(pu_col_sym); % average unit-cost of 
                                                   % rebar assembly, by 
                                                   % default
bp=b-2*rec(1);
hp=h-2*rec(2);
            
ndiam=length(RebarAvailable(:,1));
nopciones=0;
while nopciones==0
    bestArea=inf;
    maxef=1.0;
    bestMr=[];
    for i=1:ndiam % for each type of rebar
        op=i;
        
        ov=RebarAvailable(i,1); % rebar's eight-of-an-inch
        dv=RebarAvailable(i,2); % rebar's diameter
        av=(dv)^2*pi/4; % rebar's cross-sectional area
        
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
        if fdpc<2000
            sepMin=max([1.5*2.54, 1.5*dv]); % min rebar separation (cm)
        else
            sepMin=max([1.5, 1.5*dv]); % min rebar separation (in)
        end
        % There is a limit of the number of rebars that can be laid out
        % on each boundary of the cross-section, for each type of rebar
        maxVarillasSup=fix(bp/(sepMin+dv))+1;
        maxVarillasCos=fix(hp/(sepMin+dv))-1;

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            continue;
        elseif (2*maxVarillasSup)>nv
            continue;
        else
            bestCost1=nv*av*wac*height*pu_col_sym(i);
            
            for type=minVarillasSup:maxVarillasSup
                varSup=type;
                varCos=0.5*(nv-2*varSup);

                % Symmetrical design with only one type of rebar
                arregloVar1=zeros(nv,1)+i;
                [disposicion_varillado]=RebarDisposition(b,...
                     h,rec,dv,nv,varCos,varSup);
                if nv~=length(disposicion_varillado)
                    break;
                end
                % Symmetrical design with only one type of rebar
                [bestdiagram1,bestEf1,eficiencia,bestcxy1]=...
                    diagramasDisposicion(ast,b,h,E,npdiag,fdpc,nv,beta1,...
                    ov,av,disposicion_varillado,load_conditions);
                
                nvxy=[varSup varCos];
                arraySymOriginal=[varSup varSup varCos varCos];
                ab1=varSup*av;
                ab2=varCos*av;
                bestasbar1=ast;
                bestnv1=nv;
                av4_1=[ab1 ab1 ab2 ab2];
                nv4_1=arraySymOriginal;
                bestMr1=[eficiencia(1,5),eficiencia(1,7)];
                bestDisposition1=disposicion_varillado;
                
                % Symmetrical design with as many as 2 types of rebar
                [av4_2,relyEffList,bestasbar2,bestEf2,bestdiagram2,arregloVar2,...
                bestDisposition2,bestMr2,bestcxy2,bestCost2]=sym2typeRebar...
                (disposicion_varillado,op,arraySymOriginal,RebarAvailable,b,h,...
                fy,fdpc,beta1,E,load_conditions,wac,height,npdiag,ductility,...
                puSym2cols);
                
                bestnv2=nv;
                nv4_2=arraySymOriginal;
                
                % Comparison of best solutions
                if isempty(bestasbar2)==0 && isempty(bestasbar1)==0
                    
                    nopciones=nopciones+1;
                    if bestasbar2<bestasbar1
                        if bestasbar2<bestArea
                            bestdiagram=bestdiagram2;
                            bestDisposition=bestDisposition2;
                            bestArrangement=arregloVar2;
                            bestArea=bestasbar2;
                            bestCost=bestCost2;
                            
                            bestEf=bestEf2;
                            bestMr=bestMr2;
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
                            bestCost=bestCost1;
                            
                            bestEf=bestEf1;
                            bestMr=bestMr1;
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
                elseif isempty(bestasbar2)==1 && isempty(bestasbar1)==0
                    
                    if bestasbar1<bestArea && bestEf1<=maxef
                        nopciones=nopciones+1;
                        
                        bestdiagram=bestdiagram1;
                        bestDisposition=bestDisposition1;
                        bestArrangement=arregloVar1;
                        bestArea=bestasbar1;
                        bestCost=bestCost1;
                        
                        bestEf=bestEf1;
                        bestMr=bestMr1;
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
                    else
                        
                    end
                elseif isempty(bestasbar2)==1 && isempty(bestasbar1)==1
                    continue;
                end
            end
        end
    end

    if nopciones==0
        fprintf('\nThe cross-section is too small for rebar.\n');
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
        Inertia_xy_modif=[];
        break;
    else
        Mr_col=bestMr;
        if length(bestArrangement(1,:))>1
            bestArrangement=bestArrangement';
        end
        
        % Computation of cross-section's modified inertia:

        % computation of reinforcing area on each cross-section's boundary
        area_vx1=vx1Ec*av1Ec; 
        area_vx2=vx2Ec*av2Ec; 
        area_vy1=vy1Ec*av3Ec;
        area_vy2=vy2Ec*av4Ec;

        % equivalent ISR's width for each of the cross-section's boundaries
        t1_var=area_vx1/(b-2*rec(1));
        t2_var=area_vx2/(b-2*rec(1));
        t3_var=area_vy1/(h-2*rec(2));
        t4_var=area_vy2/(h-2*rec(2));

        % max load eccentricities
        [pu,imaxP]=max(abs(load_conditions(:,2)));
        pu=sign(load_conditions(imaxP,2))*pu;
        mux=max(abs(load_conditions(:,3)));
        muy=max(abs(load_conditions(:,4)));

        excentricity_x=abs(mux/pu);
        excentricity_y=abs(muy/pu);
        eccentricity_xy=[excentricity_x,excentricity_y];

        [Inertia_xy_modif,Atransf_xy,elimxy]=CrackingColumnsSym(h,b,fdpc,rec,...
            t1_var,eccentricity_xy,t3_var,pu,bestcxy,condition_cracking,E);
    end
end
% Ploteo de Diagramas de Interacción y detallado de refuerzo .........
if plotRebarDesign==1 && nopciones>0
    diagramsFinalRebarCols(load_conditions,bestdiagram,bestDisposition,...
                h,b,bestArrangement);
end