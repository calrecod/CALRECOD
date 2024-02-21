function [h,bestArea,bestEf,bestdiagram,bestArrangement,bestDisposition,...
    bestMr,bestcxy,bestCP,bestCost,bestLoad,bestCFA]=AssemRebarSymAsymIntSurf(b,...
    h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,wac,RebarAvailable,...
    load_conditions,ductility,puCostCardBuild,dataCFA,plotRebarDesign)
%-------------------------------------------------------------------------
% Syntax:
% [h,bestArea,bestEf,bestdiagram,bestArrangement,bestDisposition,...
%  bestMr,bestcxy,bestCP,bestCost,bestLoad,bestCFA]=AssemRebarSymAsymIntSurf(b,...
%  h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,wac,RebarAvailable,...
%  load_conditions,ductility,puCostCardBuild,dataCFA,plotRebarDesign)
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
%         plotRebarDesign:      is the parameters that indicates if the 
%                               rebar design results are required or not. 
%                               Options are: (1) they are required, 
%                               (2) they are not required
%
%         puCostCardBuild:      is a vector containing the parameters
%                               required for the calculation of the unit
%                               cost of a rebar design with a 
%                               "unitCostCardColsRec"
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

iter=0; maxiter=1;
noptions=0;
while noptions==0
    iter=iter+1;
    bestArea=inf;
    % --------------------------------------------------------------------
    % RP-5 (Asym-Sym4Diam)
    % --------------------------------------------------------------------

    % Optimal asymmetrical design in individual rebars (option avilable:
    % as many as four rebar diameters symmetrically distributed in number) 
    [Mr_col1,h,bestArea1,cost_elem_col1,bestdiagram1,bestnv1,Ef_sec_col1,...
    bestArrangement1,bestDisposition1,nv4,bestcxy1,bestCP1,bestLoad1,bestCFA1]=...
    AsymmSym4DiamIntSurf(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,height,...
    wac,RebarAvailable,load_conditions,ductility,puCostCardBuild,dataCFA);

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
            bestCFA=bestCFA1;
        end
    end
    % --------------------------------------------------------------------
    % RP-6 (Asym-1Diam)
    % --------------------------------------------------------------------
    
    [Mr_col2,h,bestArea2,bestCost2,bestdiagram21,bestnv2,Ef_sec_col2,...
    bestArrangement2,bestDisposition2,nv42,bestcxy2,bestCP2,bestLoad2,bestCFA2]=...
    Asym1DiamIntSurf(b,h,rec,Ac_sec_elem,npdiag,fdpc,beta1,load_conditions,...
    height,wac,RebarAvailable,ductility,puCostCardBuild,dataCFA);

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
            bestCFA=bestCFA2;
        end
    end

    % --------------------------------------------------------------------
    % RP-7 (Asym-4Diam)
    % --------------------------------------------------------------------

    [Mr_col3,h,bestArea3,bestCost3,bestdiagram31,bestnv3,Ef_sec_col3,...
    bestArrangement3,bestDisposition3,nv43,bestcxy3,bestCP3,bestLoad3,bestCFA3]=...
    Asym4DiamIntSurf(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,beta1,RebarAvailable,...
    height,wac,load_conditions,ductility,puCostCardBuild,dataCFA);

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
            bestCFA=bestCFA3;
        end
    end
    
    if (isempty(Mr_col1)==1 && isempty(Mr_col2)==1 && ...
            isempty(Mr_col3)==1)
        
        if iter <= maxiter
            if fdpc<2000 % units: kg,cm
                h=h+5;
            else       % units: lb,in
                h=h+2;
            end
            continue;
        else
            bestArea=[];
            bestEf=[];
            bestdiagram=[];
            bestArrangement=[];
            bestDisposition=[];
            bestMr=[];
            bestcxy=[];
            bestCP=[];
            bestCost=[];
            bestLoad=[];
            bestCFA=[];
            break;
        end

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