function [bestav4,relyEffList,bestArea,bestEf,bestdiagram,bestdiagram2,...
    bestArrangement,bestDisposition,bestMr,bestcxy,bestCP,bestCost]=...
    asymSym4Diam(OriginalDisposition,op,arrayOriginal,RebarAvailable,rec,...
    b,h,fy,fdpc,beta1,E,pu_asym_cols,height,wac,load_conditions,npdiag,...
    ductility)

%-------------------------------------------------------------------------
% Syntax:
% [bestav4,relyEffList,bestArea,bestEf,bestdiagram,bestdiagram2,...
%  bestArrangement,bestDisposition,bestMr,bestcxy,bestCP,bestCost]=...
%  asymSym4Diam(OriginalDisposition,op,arrayOriginal,RebarAvailable,rec,...
%  b,h,fy,fdpc,beta1,E,pu_asym_cols,height,wac,load_conditions,npdiag,...
%  ductility)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%-------------------------------------------------------------------------
% PURPOSE: To determine an optimal asymmetrical arrangement of rebars
% over a column cross-section through linear search. The distributions of
% rebars are symmetrical in quantity of rebars but asymmetrical in rebar 
% diameters. As many as four rebar diameters are allowed to be 
% simulataneously placed.
% 
% OUTPUT: bestMr:               are the final resistant bending moment for
%                               both axis directions of the optimal designed 
%                               cross-section
%
%         bestArea:             is the optimal rebar area
%
%         bestav4;              is the array containing the rebar area at
%                               each of the four cross-section boundaries
%                               corresponding to the optimal rebar design
%
%         bestCost:             is the cost of the optimal design option
%
%         bestEf:               is the critical structural efficiency 
%                               corresponding to the optimal design option
%                               against the most critical of the given load
%                               conditions
%
%         bestArrangement:      is the list of rebar type of each rebar: 
%                               size [nbars,1] (a number from 1 to N#
%
%         bestDisposition:      is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design option
%
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
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
%         pu_asym_cols:         is the average unit cost of rebar assembly
%                               for this asymmetrical rebar prototype
%
%         OriginalDisposition:  is the array containing the local rebar
%                               coordinates over the cross-section of the 
%                               original symmetrical rebar design (an array
%                               of 2 columns [xi,yi] and "n-bars" rows
%
%         arrayOriginal:        is the original symmetrical rebar
%                               arrangement from which the resulting
%                               asymmetrical rebar designs take place. The
%                               vector contains the number of rebars at
%                               each of the four cross-section's boundaries
%                               in format: 
%                               
%            [nbars-upper, nbars-lower, nbars-left, nbars-right]
%
%         op:                   is the rebar diameter index (from the rebar
%                               database) of which the original symmetrical
%                               rebar arrangment is composed
%
%         ductility:            is the parameter that indicates which
%                               ductility demand is required for the
%                               reinforcement designs. A number between 1
%                               to 3 (see documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-21
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bestDisposition=OriginalDisposition;
fc=fdpc/0.85;
if fc<2000 % units: kg,cm - Mexican NTC-17
    if ductility==1 || ductility==2 % Ductility demand level
        aminCode=0.01*b*h;
        amax=0.06*b*h;
    elseif ductility==3
        aminCode=0.01*b*h;
        amax=0.04*b*h;
    end
else % units: lb,in - ACI 318
    aminCode=0.01*b*h;
    amax=0.04*b*h;
end
bestArea=[];
Originalnv=length(OriginalDisposition);
amin=inf;
count2=0;
list_efficiencies=[];
relyEffList=[];

% Array Variations.....................................................
variationsArray=[];
for i=1:op
    for j=1:op
        for k=1:op
            for z=1:op
                %if i==j && j==k && k==z
                %elseif i==j && k==z
                %else
                    
                count2=count2+1;
                variationsArray=[variationsArray;
                                  i,j,k,z];
                typeArray=[i,j,k,z];
                % ---------------------------------------------------
                % To analyze structural efficiency of each new array:
                % ---------------------------------------------------
                ab1=arrayOriginal(1)*RebarAvailable(i,2)^2*pi/4;
                ab2=arrayOriginal(2)*RebarAvailable(j,2)^2*pi/4;
                ab3=arrayOriginal(3)*RebarAvailable(k,2)^2*pi/4;
                ab4=arrayOriginal(4)*RebarAvailable(z,2)^2*pi/4;

                ast=ab1+ab2+ab3+ab4;

                [diagrama,diagram2,pot,poc,cp,cxy,cxy2]=...
                EvalAsymDoubleDirection(npdiag,typeArray,b,h,...
                fy,fdpc,beta1,E,arrayOriginal(1),arrayOriginal(2),...
                arrayOriginal(3),arrayOriginal(4),RebarAvailable,...
                bestDisposition,rec);

                % ----------------------------------------------------
                % Analyse the resistance efficiency of the cross-
                % section, given the load conditions:
                % ----------------------------------------------------
                [maxef,eficiencia,cxy]=effRecColsDoubleDirecLS(diagrama,...
                       diagram2,load_conditions,pot,poc,cxy,cxy2);

                list_efficiencies=[list_efficiencies;
                                    maxef];

                if maxef<1.0 && ast<amin && ast>=aminCode
                    amin=ast;
                    bestArea=amin;
                    bestav4=[ab1,ab2,ab3,ab4];

                    bestEf=maxef;
                    bestdiagram=diagrama;
                    bestdiagram2=diagram2;
                    bestArrangement=zeros(Originalnv,1);
                    bestArrangement(1:arrayOriginal(1))=...
                        bestArrangement(1:arrayOriginal(1))+i;

                    bestArrangement(1+arrayOriginal(1):arrayOriginal(1)+...
                        arrayOriginal(2))=bestArrangement(1+...
                        arrayOriginal(1):arrayOriginal(1)+...
                        arrayOriginal(2))+j;

                    bestArrangement(1+arrayOriginal(1)+...
                        arrayOriginal(2):arrayOriginal(1)+...
                        arrayOriginal(2)+arrayOriginal(3))=...
                        bestArrangement(1+arrayOriginal(1)+...
                        arrayOriginal(2):arrayOriginal(1)+...
                        arrayOriginal(2)+arrayOriginal(3))+k;

                    bestArrangement(1+arrayOriginal(1)+arrayOriginal(2)+...
                        arrayOriginal(3):arrayOriginal(1)+...
                        arrayOriginal(2)+arrayOriginal(3)+...
                        arrayOriginal(4))=bestArrangement(1+...
                        arrayOriginal(1)+arrayOriginal(2)+...
                        arrayOriginal(3):arrayOriginal(1)+...
                        arrayOriginal(2)+arrayOriginal(3)+...
                        arrayOriginal(4))+z;

                    bestCost=bestArea*height*wac*pu_asym_cols;

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
if isempty(bestArea)==1 % when no feasible rebar option is found
    bestArea=[];

    bestDisposition=[];
    bestEf=[];
    bestdiagram=[];
    bestdiagram2=[];
    bestArrangement=[];
    bestav4=[];
    bestCost=[];
    bestcxy=[];
    bestCP=[];
    bestMr=[];

end