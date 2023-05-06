function [bestav4,bestnv4,relyEffList,bestArrangement,bestDisposition,...
    bestnv,bestMr,bestEf,bestcxy,bestCP,bestArea,bestdiagram,bestdiagram2,...
    bestCost]=asym1typeRebar(fdpc,nvxy,arraySymOriginal,b,h,rec,...
    RebarAvailable,op,av,npdiag,symCost,pu_asym_cols,wac,height,...
    load_conditions,ductility,beta1)

%-------------------------------------------------------------------------
% Syntax:
% [bestav4,bestnv4,relyEffList,bestArrangement,bestDisposition,...
%  bestnv,bestMr,bestEf,bestcxy,bestCP,bestArea,bestdiagram,bestdiagram2,...
%  bestCost]=asym1typeRebar(fdpc,nvxy,arraySymOriginal,b,h,rec,...
%  RebarAvailable,op,av,npdiag,symCost,pu_asym_cols,wac,height,...
%  load_conditions,ductility,beta1)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars asymmetrically 
% distributed over a rectangular column's cross-section. Only one rebar
% diameter is allowed. 
% 
% OUTPUT: bestMr:               are the final resistant bending moment for
%                               both axis directions of the optimal designed 
%                               cross-section
%
%         bestEf:               is the critical structural efficiency
%                               corresponding to the most critical load
%                               condition found in the array
%                               "load_conditions"
%
%         bestArea:             is the optimal rebar area
%
%         bestCost:             is the cost of the optimal design option
%
%         bestnv:               is the total number of rebars over the 
%                               cross-section corresponding to the optimal
%                               design option
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
%                               "Asym-1Diam".
%
%         arraySymOriginal:     is the original symmetrical rebar arrange-
%                               ment from which the resulting asymmetrical
%                               rebar designs take place. The vector con-
%                               tains the number of rebars at each of the
%                               four cross-section boundaries in format:
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
%         symCost:              is the original construction cost of the
%                               symmetrical rebar design from which the
%                               designs here generated take place. Is used
%                               to calculate the cost savings of the
%                               asymmetrical designs here generated
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-03-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
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
    fy=4200;
else       % units: lb,in - ACI 318
    aminCode=0.01*b*h;
    amax=0.04*b*h;
    
    fy=60000; % psi
end

bestArea=[];
E=fy/0.0021; % Modulus of Elasticity of reinforcing steel
amin=inf;
count=0;
list_efficiencies=[];
relyEffList=[];

% Array Variations.....................................................
variationsArray=[];
for i=0:(nvxy(1)-2)
    for j=0:(nvxy(1)-2)
        for k=0:nvxy(2)
            for z=0:nvxy(2)
                if i==j && k==z % This conditions is to exclude
                                % symmetrical designs 
                   % do nothing !
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
                    [disposition_rebar,separation_hor1,separation_hor2,...
                    separation_ver1,separation_ver2]=dispositionRebarAsymmetric...
                    (b,h,rec,nv,arrayAsym(1),arrayAsym(2),arrayAsym(3),...
                    arrayAsym(4),RebarAvailable,op,op,op,op);
                    
                    % To analyze structural efficiency of each new array
                    ast=nv*av;
                    [diagrama,diagram2,pot,poc,cp,cxy,cxy2]=...
                    EvalAsymDoubleDirection(npdiag,[op,op,op,op],b,h,...
                    fy,fdpc,beta1,E,arrayAsym(1),arrayAsym(2),arrayAsym(3),...
                    arrayAsym(4),RebarAvailable,disposition_rebar,rec);
                    
                    % Analyse the resistance efficiency of the cross-
                    % section, given the load conditions ----------------
                    [maxef,eficiencia,cxy]=effRecColsDoubleDirecLS(diagrama,...
                           diagram2,load_conditions,pot,poc,cxy,cxy2);
                    
                    list_efficiencies=[list_efficiencies;
                                        maxef];
                                    
                    if maxef<1.0 && ast<amin && ast>=aminCode
                        amin=ast;
                        bestArea=amin;
                        
                        bestDisposition=disposition_rebar;
                        bestEf=maxef;
                        bestdiagram=diagrama;
                        bestdiagram2=diagram2;
                        bestnv=nv;
                        bestArrangement=zeros(nv,1)+op;
                        bestnv4=arrayAsym;
                        bestav4=[ab1,ab2,ab3,ab4];
                        
                        bestCost=bestArea*height*wac*pu_asym_cols;
                        CostSaving=symCost-bestCost;
                        bestcxy=cxy;
                        bestCP=cp;
                        
                        bestMrx=eficiencia(1,5);
                        bestMry=eficiencia(1,7);
                        bestMr=[bestMrx bestMry];
                    end
                    if maxef<1.0
                        relyEffList=[relyEffList; maxef];
                    end
                end
            end
        end
    end
end
if isempty(bestArea)==1 % if no reliable option was found
    bestArea=[];

    bestDisposition=[];
    bestEf=[];
    bestdiagram=[];
    bestdiagram2=[];
    bestnv=[];
    bestArrangement=[];
    bestnv4=[];
    bestav4=[];

    bestCost=[];
    bestcxy=[];
    bestCP=[];

    bestMr=[];
end
