function [bestav4,relyEffList,bestArea,bestEf,bestdiagram,bestArrangement,...
    bestDisposition,bestMr,bestcxy,bestCost]=sym2typeRebar(ObarDisposition,...
    op,arraySym,RebarAvailable,b,h,fy,fdpc,beta1,E,load_conditions,wac,...
    height,npdiag,ductility,pu_sym2_cols)

%-------------------------------------------------------------------------
% Syntax:
% [bestav4,relyEffList,bestArea,bestEf,bestdiagram,bestArrangement,...
%  bestDisposition,bestMr,bestcxy,bestCost]=sym2typeRebar(ObarDisposition,...
%  op,arraySym,RebarAvailable,b,h,fy,fdpc,beta1,E,load_conditions,wac,...
%  height,npdiag,ductility,pu_sym2_cols)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design over a 
% rectangular column cross-section. As many as two rebar diameters are
% allowed. 
% 
% OUTPUT: bestMr:               are the final resistant bending moment for
%                               both axis directions of the optimal designed 
%                               cross-section
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
%                               condition
%
%         bestav4:              is the array containing the rebar area at 
%                               each of the four cross-section boundaries
%                               corresponding to the optimal rebar design
%
%         bestArrangement:      is the list of rebar type of each rebar: 
%                               size [nbars,1] (a number from 1 to 7 by 
%                               default)
%
%         bestDisposition:      is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design
%
%         bestdiagram:          is the interaction diagram data of the 
%                               optimal rebar design (considering only 
%                               positive bending moments)
%
%         bestcxy:              is a vector containing the neutral axis
%                               depth values corresponding to the
%                               most critical load condition for each of
%                               the two cross-section axis
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
%         ObarDisposition:      is the array containing the local rebar
%                               coordinate positions of the original 
%                               symmetrical rebar design from which the 
%                               symmetrical permutaitons take place  
%
%         arraySym:             is the original symmetrical rebar arrange-
%                               ment from which the resulting asymmetrical
%                               rebar designs take place. The vector
%                               arraySym the number of rebars at each of
%                               the four cross-section's boundaries in
%                               format:
%
%              [nbars-upper, nbars-lower, nbars-left, nbars-right]
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
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-21
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bestArea=[];
bestDisposition=ObarDisposition;
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

Originalnv=sum(arraySym);
amin=inf;
count2=0;
relyEffList=[];

%% Rebar permutations
variationsArray=[];
for i=1:op
    for j=1:op     
        if i~=j
            count2=count2+1;
            variationsArray=[variationsArray;
                              i,i,j,j];
            typeArray=[i,i,j,j];

            % To analyze structural efficiency of each new array
            ab1=arraySym(1)*RebarAvailable(i,2)^2*pi/4;
            ab2=arraySym(2)*RebarAvailable(i,2)^2*pi/4;
            ab3=arraySym(3)*RebarAvailable(j,2)^2*pi/4;
            ab4=arraySym(4)*RebarAvailable(j,2)^2*pi/4;

            ast=ab1+ab2+ab3+ab4;

            [maxef,diagrama,eficiencia,cp_axis,cxy]=EvaluateAsymmetric...
            (load_conditions,npdiag,typeArray,b,h,fy,fdpc,beta1,...
            E,arraySym(1),arraySym(2),arraySym(3),arraySym(4),...
            RebarAvailable,ObarDisposition);
        
            if maxef<1.0 && ast<amin && ast>=aminCode
                amin=ast;
                bestArea=amin;
                bestav4=[ab1,ab2,ab3,ab4];

                bestEf=maxef;
                bestdiagram=diagrama;
                bestArrangement=zeros(1,Originalnv);
                
                bestArrangement(1:2*arraySym(1))=...
                    bestArrangement(1:2*arraySym(1))+i;

                bestArrangement(1+2*arraySym(1):2*arraySym(1)+2*arraySym(3))=...
                    bestArrangement(1+2*arraySym(1):2*arraySym(1)+2*...
                    arraySym(3))+j;

                bestCost=bestArea*height*wac*pu_sym2_cols;

                bestcxy=cxy;

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
if isempty(bestArea)==1

    bestDisposition=[];
    bestEf=[];
    bestdiagram=[];
    bestnv=[];
    bestArrangement=[];
    bestnv4=[];
    bestav4=[];

    bestCost=[];
    bestcxy=[];

    bestMr=[];
end