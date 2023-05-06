function [Mr_col,h,Inertia_xy_modif,bestArea,bestCost,bestdiagram,...
    bestdiagram2,bestnv,bestEf,bestArrangement,bestDisposition,nv4,...
    bestcxy,bestCP]=AsymmSym4Diam(b,h,rec,act,E,npdiag,fdpc,beta1,...
    height,wac,RebarAvailable,load_conditions,pu_asym_cols,...
    condition_cracking,ductility)

%-------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,Inertia_xy_modif,bestArea,bestCost,bestdiagram,...
%  bestdiagram2,bestnv,bestEf,bestArrangement,bestDisposition,nv4,...
%  bestcxy,bestCP]=AsymmSym4Diam(b,h,rec,act,E,npdiag,fdpc,beta1,...
%  height,wac,RebarAvailable,load_conditions,pu_asym_cols,...
%  condition_cracking,ductility)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal asymmetrical arrangement of rebars 
% over a column cross-section through linear search. Every single feasible
% symmetrical reinforcement option is assessed, from which the asymmetrical
% design then takes place. As many as four rebar diameters are allowed to 
% be placed simultaneously asymmetrically, although symmetrically 
% distributed in number for each cross-section boundary.
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
%         bestnv:               is the total number of rebars over the 
%                               cross-section corresponding to the optimal
%                               design option
%
%         bestEf:               is the critical structural efficiency 
%                               corresponding to the optimal rebar option
%                               against the most critical of the given load
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
%                               optimal rebar design (considering only 
%                               positive bending moments)
%
%         bestdiagram2:         is the interaction diagram data of the 
%                               optimal rebar design (considering only 
%                               negative bending moments)
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
%         pu_asym_cols:         is the average unit construction cost for
%                               this rebar prototype
%
%         condition_cracking:   parameter that indicates which cross-section
%                               cracking mechanism will be consider, either 
%                               Cracked or Non-cracked. If the condition 
%                               Non-cracked is set, then the cracking 
%                               mechanism will be neglected by all means
%
%         ductility:            is a parameter that indicates which
%                               ductility demand level is required for the
%                               reinforcement designs. A number between 1
%                               to 3 (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-21
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
fc=fdpc/0.85;
fy=E*0.0021; % yield stress of reinforcing steel
bp=b-2*rec(1);
hp=h-2*rec(2);

ndiam=length(RebarAvailable(:,1));
noptions=0;
while noptions==0

    bestArea=inf;
    for i=1:ndiam % for each type of rebar
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
                
                % Asymmetrical design with as many as 4 types of rebar
                [av4_2,relyEffList,bestasbar2,bestEf2,bestdiagram21,...
                bestdiagram22,arregloVar2,bestDisposition2,bestMr2,...
                bestcxy2,bestCP2,bestCost2]=asymSym4Diam...
                (disposicion_varillado,op,arraySymOriginal,RebarAvailable,...
                rec,b,h,fy,fdpc,beta1,E,pu_asym_cols,height,wac,...
                load_conditions,npdiag,ductility);
                
                bestnv2=nv;
                nv4_2=arraySymOriginal;
                
                if isempty(bestasbar2)==0
                    noptions=noptions+1;
                    if bestasbar2<bestArea
                        bestdiagram=bestdiagram21;
                        bestdiagram2=bestdiagram22;
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
                    
                elseif isempty(bestasbar2)==1
                    continue;
                end
                   
            end
        end
    end
    
    if noptions==0
        fprintf('\nDimensions of column are too small.\n');
        fprintf('The cross-section dimensions should be increased.\n');

        bestdiagram=[];
        bestdiagram2=[];
        bestDisposition=[];
        bestArrangement=[];
        bestArea=[];
        bestCost=[];
        bestEf=[];
        bestMr=[];
        bestnv=[];
        nv4=[];
        av4=[];
        bestcxy=[];
        bestCP=[];
        Inertia_xy_modif=[];
        Mr_col=bestMr;
        
        break;
    else
        
        Mr_col=bestMr;
        %% Computation of cross-section's modified inertia:

        % computation of reinforcing area on each cross-section's 
        % boundaries
        area_vx1=vx1Ec*av1Ec; 
        area_vx2=vx2Ec*av2Ec; 
        area_vy1=vy1Ec*av3Ec;
        area_vy2=vy2Ec*av4Ec;

        % equivalent ISR's width for each of the four cross-section's
        % boundaries
        t1_var=area_vx1/(b-2*rec(1));
        t2_var=area_vx2/(b-2*rec(1));
        t3_var=area_vy1/(h-2*rec(2));
        t4_var=area_vy2/(h-2*rec(2));

        % Max load eccentricities for each cross-section's axis
        pu=load_conditions(1,2);
        mux=max(abs(load_conditions(:,3)));
        muy=max(abs(load_conditions(:,4)));

        excentricity_x=abs(mux/pu);
        excentricity_y=abs(muy/pu);
        eccentricity_xy=[excentricity_x,excentricity_y];

        % Cracked/Non-cracked transoformed section mechanism
        [Inertia_xy_modif,Atransfxy]=CrackingColumnsAsym(h,b,fdpc,rec,...
            eccentricity_xy,t1_var,t2_var,t3_var,t4_var,pu,bestcxy,...
            condition_cracking,bestCP);

    end
end