function [Mr_col,h,bestArea,lowestCost,ovMostEc,nvEc,maxEfEc,...
    bestArrangement,bestDisposition,nvxy,bestdiagram,bestLoad]=...
    optRebarColsSym2PackIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,...
    RebarAvailable,wac,height,load_conditions,puCostCardBuild,...
    plotRebarDesign)

%------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,bestArea,lowestCost,ovMostEc,nvEc,maxEfEc,...
% bestArrangement,bestDisposition,nvxy,bestdiagram,bestLoad]=...
% optRebarColsSym2PackIntSurf(b,h,rec,act,E,npdiag,fdpc,beta1,...
% RebarAvailable,wac,height,load_conditions,puCostCardBuild,...
% plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design composed of
%          only one type of rebar in packages of two rebars.
% 
% Note: The structural efficiency for each rebar design is determined by
% rotating the cross-section according to each given load combination, so
% that an interaction diagram is computed for each load combination's
% direction of action.
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
%         lowestCost:           is the cost of the optimal design option
%
%         ovMostEc:             is the type of rebar corresponding to the
%                               optimal rebar design option
%
%         nvEc:                 is the total number of rebars over the 
%                               cross-section corresponding to the optimal
%                               design option
%
%         maxEfEc:              is the critical structural efficiency 
%                               corresponding to the optimal most economic
%                               design option maxEfEc<1.0
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
%         nvxy:                 number of rebars placed horizontally and
%                               vertically [nbar-x,nbars-y] 
%
%         bestLoad:             is the resultant most critical load 
%                               combination, in format: [n-load, Pu, Mu],
%                               where Mu = sqrt( Mux^2 + Muy^2 )
%
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
%
%         act:                  optima ISR reinforcement area
%
%         sepMin:               min separation of rebars constraint
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
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
pucb=puCostCardBuild;
pu_col_sym=unitCostCardColsRec(pucb(1),pucb(2),pucb(3),pucb(4),pucb(5),...
                               pucb(6),pucb(7));
fy=E*0.0021;

bp=b-2*rec(1);
hp=h-2*rec(2);

ntrebar=length(RebarAvailable(:,1));
noptions=0;
while noptions==0
    
    bestArea=inf;
    maxEfEc=1;
    for i=1:ntrebar
        
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
        % Rebar separation:
        if fdpc<2000  % units (kg,cm)
            sepMin=max([1.5*2.54, 3/2*dv]);
        else          % units (lb,in)
            sepMin=max([1.5, 3/2*dv]);
        end
        % There is a limit of the number of rebars that can be laid out
        % on each boundary of the cross-section, for each type of rebar
        maxVarillasSup=2*(fix((bp)/(sepMin+2*dv)))+2;
        maxVarillasCos=2*(fix((hp)/(sepMin+2*dv)));
        if 2*maxVarillasSup>nv
            maxVarillasSup=nv/2;
        end

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        if 2*maxVarillasCos>(nv-4)
            maxVarillasCos=(nv-4)/2;
        end
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            continue;
        else
            
            costo=nv*av*height*wac*pu_col_sym;
            for nvsup=minVarillasSup:maxVarillasSup
                varSup=nvsup;
                varCos=0.5*(nv-2*varSup);
                arregloVar=zeros(nv,1)+i;
                [dispRebar]=RebarDisposition2packSym(b,...
                                        h,rec,dv,nv,varCos,varSup);
                
                [eficiencia,iloadmax,bestLoad1,gamma,diagrama,...
                rotdispositionRebar,rotsection,cxy,cp]=multiDiagAxisColRec...
                (b,h,load_conditions,[i,i,i,i],npdiag,fy,fdpc,beta1,E,...
                varSup,varSup,varCos,varCos,RebarAvailable,dispRebar);
                
                maxef=eficiencia(iloadmax,5);
                
                if ast<=bestArea && maxef<maxEfEc
                    nvEc=nv;
                    lowestCost=costo;
                    ovMostEc=ov;
                    bestLoad=bestLoad1;
                    Mr_col=eficiencia(iloadmax,4);
                    maxEfEc=maxef;
                    nvxy=[varSup,varCos];
                    bestDisposition=dispRebar;
                    bestdiagram=diagrama;
                    bestArrangement=arregloVar;
                    bestArea=ast;
                    best_c=cxy;
                    
                    noptions=noptions+1;
                end
            end
        end
    end
    
    if noptions==0
        fprintf('\nThe dimensions of the column are not fit for rebar.\n');
        fprintf('The dimensions should be increased.\n');

        % h=h+5;
        
        nvEc=nv;
        lowestCost=[];
        ovMostEc=[];

        maxEfEc=[];
        bestDisposition=[];
        bestdiagram=[];
        bestLoad=[];
        bestArrangement=[];
        bestArea=[];
        nvxy=[];
        Mr_col=[];
        
        break;
    end
end
if plotRebarDesign==1 && noptions>0
    section=[0.5*b 0.5*h;
             -0.5*b 0.5*h;
             -0.5*b -0.5*h;
              0.5*b -0.5*h;
              0.5*b 0.5*h];
      
    PlotRotRecSecRebarCols(section,bestLoad,bestdiagram,...
                    bestDisposition,b,h,bestArrangement);
end