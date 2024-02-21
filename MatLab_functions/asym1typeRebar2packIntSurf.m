function [bestav4,bestnv4,bestArrangement,bestDisposition,...
    bestnv,bestMr,bestEf,bestcxy,bestCP,bestArea,bestdiagram,...
    bestCost,maxLoad,bestCFA]=asym1typeRebar2packIntSurf(fdpc,fy,nvxy,...
    arraySymOriginal,b,h,rec,RebarAvailable,op,av,npdiag,...
    height,wac,load_conditions,ductility,beta1,puCostCardBuild,dataCFA)

%-------------------------------------------------------------------------
% Syntax:
% [bestav4,bestnv4,bestArrangement,bestDisposition,...
% bestnv,bestMr,bestEf,bestcxy,bestCP,bestArea,bestdiagram,...
% bestCost,maxLoad,bestCFA]=asym1typeRebar2packIntSurf(fdpc,fy,nvxy,...
% arraySymOriginal,b,h,rec,RebarAvailable,op,av,npdiag,symCost,...
% height,wac,load_conditions,ductility,beta1,puCostCardBuild,dataCFA)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars asymmetrically 
% distributed over a rectangular column's cross-section in packages of two
% rebars. Only one rebar diameter is allowed for each design option.
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
%                               against the most critical of the load
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
%         maxLoad:              is the resultant most critical load 
%                               combination, in format: [n-load, Pu, Mu],
%                               where Mu = sqrt( Mux^2 + Muy^2 )
%
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
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
%                               this rebar prototype (corresponding to the
%                               "Asym-2pack-1Diam".
%
%         arraySymOriginal:     is the original symmetrical rebar arrange-
%                               ment with packages of two rebars, from 
%                               which the resulting asymmetrical rebar
%                               designs take place. The vector contains
%                               the number of rebars at each of the four
%                               cross-section's boundaries in format:
%
%              [nbars-upper, nbars-lower, nbars-left, nbars-right]
%
%         nvxy:                 is a vector containing the number of rebars
%                               in the upper or lower boundary of the
%                               cross-section and in the left or right
%                               boundary of the cross-section, in format:
%                               [nbars-upper, nbars-lower]
%
%         op:                   is the rebar diameter index (from the rebar
%                               database table - a number between 1 to 7)
%                               of which the rebar design is composed
%
%         ductility:            is a parameter that indicates which
%                               ductility demand is reuired for the
%                               reinforcement designs. A number between 1
%                               to 3 (see Documentation).
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
fc=fdpc/0.85;
if fc<2000 % units: kg,cm - Mexican NTC-17
    if ductility==1 || ductility==2 % ductility demand level
        aminCode=0.01*b*h;
        amax=0.06*b*h;
    elseif ductility==3
        aminCode=0.01*b*h;
        amax=0.04*b*h;

    end
else       % units: lb,in - ACI 318
    aminCode=0.01*b*h;
    amax=0.04*b*h;
end
bestArea=[];
E=fy/0.0021;
amin=inf;
count=0;

% Array Combinations
variationsArray=[];
for i=0:(nvxy(1)-2)
    for j=0:(nvxy(1)-2)
        for k=0:nvxy(2)
            for z=0:nvxy(2)
                if i==j && k==z % This condition is to exclude
                                 % symmetrical designs
                    variationsArray=variationsArray;
                else
                    count=count+1;
                    variationsArray=[variationsArray;
                                    i,j,k,z];
                    arrayAsym=arraySymOriginal-[i,j,k,z];
                    
                    ab1=arrayAsym(1)*av;
                    ab2=arrayAsym(2)*av;
                    ab3=arrayAsym(3)*av;
                    ab4=arrayAsym(4)*av;
                    
                    nv=sum(arrayAsym);
                    
                    % To re-distribute rebars over the cross-sections 
                    [dispositionRebar]=RebarDisposition2packAsym(b,...
                    h,rec,nv,arrayAsym(1),arrayAsym(2),arrayAsym(3),...
                    arrayAsym(4),RebarAvailable,op,op,op,op);
                    
                    % To analyze structural efficiency of each new array
                    ast=nv*av;
                    % ----------------------------------------------------
                    % Analyse the resistance efficiency of the cross-
                    % section, given the load conditions:
                    % ----------------------------------------------------
                    [eficiencia,iloadmax,trans_load_condition,gamma,...
                    diagrama,rotdispositionRebar,rotsection,cxy,cp]=...
                    multiDiagAxisColRec(b,h,load_conditions,[op,op,op,op],...
                    npdiag,fy,fdpc,beta1,E,arrayAsym(1),arrayAsym(2),...
                    arrayAsym(3),arrayAsym(4),RebarAvailable,dispositionRebar);

                    wub=dataCFA(3);
                    wnd=dataCFA(4);
                    [BS,CFA]=BuildabilityScoreRebarCols([op,op,op,op],...
                        arrayAsym,wub,wnd);

                    pccb=puCostCardBuild;
                    pu_asym_cols=unitCostCardColsRec(pccb(1),pccb(2),pccb(3),...
                                 pccb(4)*CFA,pccb(5),pccb(6),pccb(7));
                             
                    maxef=eficiencia(iloadmax,5);
                    
                    if maxef<1.0 && ast<amin && ast>=aminCode && ...
                            dataCFA(1)<=CFA && CFA<=dataCFA(2)
                        amin=ast;
                        bestArea=amin;
                        
                        bestDisposition=dispositionRebar;
                        bestEf=maxef;
                        bestdiagram=diagrama;
                        maxLoad=trans_load_condition;
                    
                        bestnv=nv;
                        bestArrangement=zeros(1,nv)+op;
                        bestnv4=arrayAsym;
                        bestav4=[ab1,ab2,ab3,ab4];
                        
                        bestCost=bestArea*height*wac*pu_asym_cols;
                        
                        bestcxy=cxy;
                        bestCP=cp;
                        bestCFA=CFA;
                        bestMr=eficiencia(1,4);
                    end
                end
            end
        end
    end
end
if isempty(bestArea)==1 % if no feasible option was found
    bestArea=[];

    bestDisposition=[];
    bestEf=[];
    bestdiagram=[];
    bestnv=[];
    bestArrangement=[];
    bestnv4=[];
    bestav4=[];
    bestCFA=[];
    bestCost=[];
    bestcxy=[];
    bestCP=[];

    bestMr=[];
    maxLoad=[];
end