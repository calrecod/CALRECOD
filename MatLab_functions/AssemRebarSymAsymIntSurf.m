function [h,bestArea,bestEf,bestdiagram,bestArrangement,bestDisposition,...
    bestMr,bestcxy,bestCP,bestCost,bestLoad]=AssemRebarSymAsymIntSurf(b,...
    h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,wac,pu_sym_cols,...
    RebarAvailable,load_conditions,ductility,plotRebarDesign)

%-------------------------------------------------------------------------
% Syntax:
% [h,bestArea,bestEf,bestdiagram,bestArrangement,bestDisposition,...
%  bestMr,bestcxy,bestCP,bestCost,bestLoad]=AssemRebarSymAsymIntSurf(b,...
%  h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,wac,pu_sym_cols,...
%  RebarAvailable,load_conditions,ductility,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars, either 
% symmetrical or asymmetrical over a column cross-section in individual
% rebars.
% 
% NOTE: The structural efficiency of each rebar design is determined by
% rotating the cross-section according to each given load combination so
% that their corresponding interaction diagram is computed with respect to
% the action axis of each load condition.
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
%         bestLoad:             is the resultant most critical load 
%                               combination, in format: [n-load, Pu, Mu],
%                               where Mu = sqrt( Mux^2 + Muy^2 )
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
%                               and construction unit cost for the most
%                               symmetrical design prototype: format by
%                               default:
%    -----------------------------------------------------------------
%    pu_col=[PU{#4}, PU{#5}, PU{#6}, PU{#8}, PU{#9}, PU{#10}, ...]
%    -----------------------------------------------------------------
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
    [Mr_col1,h,bestArea1,cost_elem_col1,bestdiagram1,bestnv1,Ef_sec_col1,...
    bestArrangement1,bestDisposition1,nv4,bestcxy1,bestCP1,bestLoad1]=...
    AsymmSym4DiamIntSurf(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,...
    wac,RebarAvailable,load_conditions,pu_asym_cols,ductility);

    if isempty(bestArea1)==0
        if bestArea1<bestArea
            bestArea=bestArea1;
            bestEf=Ef_sec_col1;
            bestdiagram=bestdiagram1;
            bestArrangement=bestArrangement1;
            bestDisposition=bestDisposition1;
            bestMr=Mr_col1;
            bestcxy=bestcxy1;
            bestCP=bestCP1;
            bestCost=cost_elem_col1;
            bestLoad=bestLoad1;
        end
    end
    % --------------------------------------------------------------------
    % RP-6 (Asym-1Diam)
    % --------------------------------------------------------------------
    
    pu_asym_cols=1/0.8*sum(pu_sym_cols)/length(pu_sym_cols);
    
    [Mr_col2,h,bestArea2,bestCost2,bestdiagram21,bestnv2,Ef_sec_col2,...
    bestArrangement2,bestDisposition2,nv42,bestcxy2,bestCP2,bestLoad2]=...
    Asym1DiamIntSurf(b,h,rec,Ac_sec_elem,npdiag,fdpc,beta1,load_conditions,...
    pu_asym_cols,pu_sym_cols,height,wac,RebarAvailable,ductility);

    if isempty(bestArea2)==0
        if bestArea2<bestArea
            bestArea=bestArea2;
            bestEf=Ef_sec_col2;
            bestdiagram=bestdiagram21;
            bestArrangement=bestArrangement2;
            bestDisposition=bestDisposition2;
            bestMr=Mr_col2;
            bestcxy=bestcxy2;
            bestCP=bestCP2;
            bestCost=bestCost2;
            bestLoad=bestLoad2;
        end
    end

    % --------------------------------------------------------------------
    % RP-7 (Asym-4Diam)
    % --------------------------------------------------------------------

    pu_asym_cols=1/0.7*sum(pu_sym_cols)/length(pu_sym_cols);

    [Mr_col3,h,bestArea3,bestCost3,bestdiagram31,bestnv3,Ef_sec_col3,...
    bestArrangement3,bestDisposition3,nv43,bestcxy3,bestCP3,bestLoad3]=...
    Asym4DiamIntSurf(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,RebarAvailable,...
    height,wac,load_conditions,pu_asym_cols,ductility);

    if isempty(bestArea3)==0
        if bestArea3<bestArea
            bestArea=bestArea3;
            bestEf=Ef_sec_col3;
            bestdiagram=bestdiagram31;
            bestArrangement=bestArrangement3;
            bestDisposition=bestDisposition3;
            bestMr=Mr_col3;
            bestcxy=bestcxy3;
            bestCP=bestCP3;
            bestCost=bestCost3;
            bestLoad=bestLoad3;
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

    section=[0.5*b 0.5*h;
         -0.5*b 0.5*h;
         -0.5*b -0.5*h;
          0.5*b -0.5*h;
          0.5*b 0.5*h];
      
    PlotRotRecSecRebarCols(section,bestLoad,bestdiagram,...
                    bestDisposition,b,h,bestArrangement);
end