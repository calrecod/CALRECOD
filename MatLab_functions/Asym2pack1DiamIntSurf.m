function [bestMr,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
    bestArrangement,bestDisposition,nv4,bestcxy,bestCP,bestLoad]=...
    Asym2pack1DiamIntSurf(b,h,rec,act,npdiag,fdpc,fy,beta1,load_conditions,...
    RebarAvailable,height,wac,pu_asym_cols,pu_col_sym,ductility)

%-------------------------------------------------------------------------
% Syntax:
% [bestMr,h,bestArea,bestCost,bestdiagram,bestnv,bestEf,...
% bestArrangement,bestDisposition,nv4,bestcxy,bestCP,bestLoad]=...
% Asym2pack1DiamIntSurf(b,h,rec,act,npdiag,fdpc,fy,beta1,load_conditions,...
% RebarAvailable,height,wac,pu_asym_cols,pu_col_sym,ductility)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars asymmetrically 
% distributed over a column cross-section through linear search in packages
% of two rebars. Only one rebar diameter is allowed. 
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
%         nv4:                  is the total number of rebars at each  
%                               boundary of the column cross-section 
%                               corresponding to the optimal design option
%
%         bestEf:               is the critical structural efficiency 
%                               corresponding to the optimal design option
%                               against the given load condition
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
%         
%         bestdiagram:          is the interaction diagram data of the 
%                               optimal rebar design
%
%         bestcxy:              is a vector containing the neutral axis 
%                               depth values corresponding to the most 
%                               critical load condition for each of the 
%                               two cross-section axis
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
%         act:                  optima ISR reinforcement area (from the
%                               function "isrColumns"
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
%         pu_asym_cols:         is the average construction unit cost for 
%                               this rebar prototype
%
%         pu_sym_cols:          is the database of reinforcement assembly
%                               and construction unit cost for the most
%                               symmetrical design prototype: format by
%                               default:
%    -----------------------------------------------------------------
%    pu_col=[PU{#4}, PU{#5}, PU{#6}, PU{#8}, PU{#9}, PU{#12}, ...]
%    -----------------------------------------------------------------
%
%         ductility:            is a parameter that indicates which 
%                               ductility demand is required for the 
%                               reinforcement designs. A number between 
%                               1 to 3 (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-21
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
fc=fdpc/0.85;
bp=b-2*rec(1);
hp=h-2*rec(2);

ndiam=length(RebarAvailable(:,1));
noptions=0;
while noptions==0
    bestArea=inf;
    for i=1:ndiam % for each type of rebar
        op=i;
        
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
        maxVarillasSup=2*(fix((bp)/(sepMin+2*dv)))+2;
        maxVarillasCos=2*(fix((hp)/(sepMin+2*dv)));

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            continue;
        elseif (2*maxVarillasSup<nv)
            continue;
        else
            costSym=nv*av*height*wac*pu_col_sym(i);
            
            for type=minVarillasSup:maxVarillasSup
                varSup=type;
                varCos=0.5*(nv-2*varSup);

                [disposicion_varillado]=RebarDisposition2packSym(b,...
                                        h,rec,dv,nv,varCos,varSup);
                if nv~=length(disposicion_varillado)
                    break;
                end
                nvxy=[varSup varCos];
                arraySymOriginal=[varSup varSup varCos varCos];
                
                % Asymmetrical design with only one rebar diameter
                % in packs of two
                [av4_1,nv4_1,arregloVar1,bestDisposition1,bestnv1,bestMr1,...
                bestEf1,bestcxy1,bestCP1,bestasbar1,bestdiagram1,bestCost1,...
                maxLoad]=asym1typeRebar2packIntSurf(fdpc,fy,nvxy,...
                arraySymOriginal,b,h,rec,RebarAvailable,op,av,npdiag,...
                costSym,pu_asym_cols,height,wac,load_conditions,ductility,...
                beta1);
            
                if bestasbar1<bestArea 
                    noptions=noptions+1;
                    
                    bestdiagram=bestdiagram1;
                    bestLoad=maxLoad;
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
                    bestCP=bestCP1;
                    
                    vx1Ec=nv4(1);
                    vx2Ec=nv4(2); 
                    vy1Ec=nv4(3);
                    vy2Ec=nv4(4);

                    av1Ec=av4(1); 
                    av2Ec=av4(2); 
                    av3Ec=av4(3);
                    av4Ec=av4(4);
                
                elseif isempty(bestasbar1)==1
                    continue;
                end
                   
            end
        end
    end

    if noptions==0
        fprintf('\nThe columns cross-section dimensions are too small\n');
        fprintf('rebar for this rebar prototype Asym-2pack-1Diam.\n');

        bestdiagram=[];
        bestDisposition=[];
        bestArrangement=[];
        bestArea=[];
        bestCost=[];
        bestEf=[];
        bestMr=[];
        bestLoad=[];
        bestnv=[];
        nv4=[];
        bestcxy=[];
        bestCP=[];
        break;
    end
end