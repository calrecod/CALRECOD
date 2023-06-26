function [Inertia_xy_modif,h,bestArea,bestEf,bestdiagram,bestdiagram2,...
    bestArrangement,bestDisposition,bestMr,bestcxy,bestCP,bestCost]=...
    AssembledRebarSymAsym2(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,...
    height,wac,pu_sym_cols,RebarAvailable,load_conditions,...
    condition_cracking,ductility,plotRebarDesign)

%-------------------------------------------------------------------------
% Syntax:
% [Inertia_xy_modif,h,bestArea,bestEf,bestdiagram,bestdiagram2,...
%  bestArrangement,bestDisposition,bestMr,bestcxy,bestCP,bestCost]=...
%  AssembledRebarSymAsym2(b,h,rec,Ac_sec_elem,E,npuntos,fdpc,beta1,...
%  height,wac,pu_sym_cols,RebarAvailable,load_conditions,...
%  condition_cracking,ductility,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars, either 
% symmetrical or asymmetrical over a column cross-section in individual
% rebars.
% 
% NOTE: The structural efficiency of each rebar design is determined with
% the Inverse Load method (Bresler's formula) and the Contour Load method.
% Thus, only one interaction diagram is computed for the whole given set of
% load conditions, for each rebar design.
%
% OUTPUT: bestMr:               are the final resistant bending moment for
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
%                               corresponding to the optimal design option 
%                               bestEf<1.0
%
%         bestArrangement:      is the list of rebar type indices of each 
%                               rebar: size [nbars,1] (a number from 1 to 
%                               7 by default)
%
%         bestDisposition:      is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design option
%
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
%
%         Ac_sec_elem:          optima ISR reinforcement area
%
%         sepMin:               min separation of rebars constraint
%
%         E:                    Modulus of Elasticity of the reinforcement
%                               steel
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
%         pu_sym_cols:          is the database of reinforcement assembly
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
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-03-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
fc=fdpc/0.85;
iter=0;
noptions=0;
while noptions==0
    iter=iter+1;
    bestArea=inf;
    % --------------------------------------------------------------------
    % RP-5 (Asym-Sym4Diam)
    % --------------------------------------------------------------------
    pu_asym_cols=1/0.8*sum(pu_sym_cols)/length(pu_sym_cols);

    % Optimal asymmetrical design in individual rebars (option avilable:
    % as many as four rebar diameters symmetrically distributed in number) 
    [Mr_col1,h,Inertia_xy_modif1,bestArea1,cost_elem_col1,bestdiagram1,...
    bestdiagram12,nv2,Ef_sec_col1,bestArrangement1,bestDisposition1,nv4,...
    bestcxy1,bestCP1]=AsymmSym4Diam(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,...
    beta1,height,wac,RebarAvailable,load_conditions,pu_asym_cols,...
    condition_cracking,ductility);

    if isempty(bestArea1)==0
        if bestArea1<bestArea
            bestArea=bestArea1;
            bestEf=Ef_sec_col1;
            bestdiagram=bestdiagram1;
            bestdiagram2=bestdiagram12;
            bestArrangement=bestArrangement1;
            bestDisposition=bestDisposition1;
            bestMr=Mr_col1;
            bestcxy=bestcxy1;
            bestCP=bestCP1;
            bestCost=cost_elem_col1;
            Inertia_xy_modif=Inertia_xy_modif1;
        end
    end
    % --------------------------------------------------------------------
    % RP-6 (Asym-1Diam)
    % --------------------------------------------------------------------
    
    pu_asym_cols=1/0.8*sum(pu_sym_cols)/length(pu_sym_cols);
    
    [Mr_col2,h,Inertia_xy_modif2,bestArea2,bestCost2,bestdiagram21,...
    bestdiagram22,bestnv2,Ef_sec_col2,bestArrangement2,bestDisposition2,...
    nv42,bestcxy2,bestCP2]=Asym1Diam(b,h,rec,Ac_sec_elem,npdiag,fdpc,beta1,...
    load_conditions,pu_asym_cols,height,wac,RebarAvailable,...
    condition_cracking,ductility);

    if isempty(bestArea2)==0
        if bestArea2<bestArea
            bestArea=bestArea2;
            bestEf=Ef_sec_col2;
            bestdiagram=bestdiagram21;
            bestdiagram2=bestdiagram22;
            bestArrangement=bestArrangement2;
            bestDisposition=bestDisposition2;
            bestMr=Mr_col2;
            bestcxy=bestcxy2;
            bestCP=bestCP2;
            bestCost=bestCost2;
            Inertia_xy_modif=Inertia_xy_modif2;
        end
    end

    % --------------------------------------------------------------------
    % RP-7 (Asym-4Diam)
    % --------------------------------------------------------------------

    pu_asym_cols=1/0.7*sum(pu_sym_cols)/length(pu_sym_cols);

    [Mr_col3,h,Inertia_xy_modif3,bestArea3,bestCost3,bestdiagram31,...
    bestdiagram32,bestnv3,Ef_sec_col3,bestArrangement3,bestDisposition3,...
    nv43,bestcxy3,bestCP3]=Asym4Diam(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,...
    beta1,RebarAvailable,height,wac,load_conditions,pu_asym_cols,...
    condition_cracking,ductility);

    if isempty(bestArea3)==0
        if bestArea3<bestArea
            bestArea=bestArea3;
            bestEf=Ef_sec_col3;
            bestdiagram=bestdiagram31;
            bestdiagram2=bestdiagram32;
            bestArrangement=bestArrangement3;
            bestDisposition=bestDisposition3;
            bestMr=Mr_col3;
            bestcxy=bestcxy3;
            bestCP=bestCP3;
            bestCost=bestCost3;
            Inertia_xy_modif=Inertia_xy_modif3;
        end
    end
    
    if (isempty(Mr_col1)==1 && isempty(Mr_col2)==1 && ...
            isempty(Mr_col3)==1)
        if fc<2000 % units: kg,cm
            h=h+5;
        else       % units: lb,in
            h=h+2;
        end
        continue;

    else
        break;
    end
end

if plotRebarDesign==1 && (isempty(bestArea1)==0  || ...
    isempty(bestArea2)==0 || isempty(bestArea3)==0)

    % Plot the interaction diagram for better assessment 
    diagDoubleDirecAsymRebarCols(load_conditions,bestdiagram,...
        bestdiagram2,bestDisposition,h,b,bestArrangement);
end