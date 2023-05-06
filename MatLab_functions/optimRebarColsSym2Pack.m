function [Mr_col,h,Inertia_xy_modif,bestArea,lowestCost,ovMostEc,nvEc,...
    maxEfEc,bestArrangement,best_disposicion,nvxy,bestdiagrama]=optimRebarColsSym2Pack...
    (b,h,rec,act,E,npuntos,fdpc,beta1,pu_col_sym,RebarAvailable,wac,height,...
    load_conditions,condition_cracking,plotRebarDesign)

%------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,Inertia_xy_modif,bestArea,lowestCost,ovMostEc,nvEc,...
%  maxEfEc,bestArrangement,best_disposicion,nvxy]=optimRebarColsSym2Pack...
%  (b,h,rec,act,E,npuntos,fdpc,beta1,pu_col_sym,RebarAvailable,wac,height,...
%  load_conditions,condition_cracking,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design composed of
%          only one type of rebar in packages of two rebars.
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
%         best_disposicion:     is an array containing the local coordinates 
%                               of position of each rebar over the cross-
%                               section corresponding to the optimal rebar
%                               design option
%
%         nvxy:                 number of rebars placed horizontally and
%                               vertically [nbar-x,nbars-y] 
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
%         npntos:               number of points to compute for the 
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
%         pu_col:               is the database of reinforcement assembly
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
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

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
            
            costo=nv*av*height*wac*pu_col_sym(i);
            for nvsup=minVarillasSup:maxVarillasSup
                varSup=nvsup;
                varCos=0.5*(nv-2*varSup);
                arregloVar=zeros(nv,1)+i;
                [disposicion_varillado]=RebarDisposition2packSym(b,...
                                        h,rec,dv,nv,varCos,varSup);
                
                [diagrama,maxef,eficiencia,cxy]=diagramasDisposicion(ast,...
                b,h,E,npuntos,fdpc,nv,beta1,ov,av,disposicion_varillado,...
                load_conditions);

                if ast<=bestArea && maxef<maxEfEc
                    nvEc=nv;
                    lowestCost=costo;
                    ovMostEc=ov;
                    bestMrx=eficiencia(1,5);
                    bestMry=eficiencia(1,7);
                    maxEfEc=maxef;
                    vxEc=varSup;
                    vyEc=varCos;
                    best_disposicion=disposicion_varillado;
                    bestdiagrama=diagrama;
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
        fprintf('The height dimension will be increased.\n');

        % h=h+5;
        
        nvEc=nv;
        lowestCost=[];
        ovMostEc=[];

        maxEfEc=[];
        best_disposicion=[];
        bestdiagrama=[];
        bestArrangement=[];
        bestArea=[];
        best_c=[];
        Mr_col=[];
        Inertia_xy_modif=[];
        
        noptions=1;
    else
        avEc=av;
        
        Mr_col(1)=bestMrx;
        Mr_col(2)=bestMry;

        %%% MODIFIED INERTIA MOMENTUM THROUGH THE TRANSFORMED CRACKED
        %%% CROSS-SECTION .............................................

        %%%% Calculation of steel area on each cross-section's boundary
        area_vx=vxEc*avEc; 
        area_vy=vyEc*avEc;
        nvxy=[vxEc,vyEc];

        t1_var=area_vx/(b-2*rec(1));
        t2_var=area_vy/(h-2*rec(2));

        % Load eccentricities
        [pu,imaxP]=max(abs(load_conditions(:,2)));
        pu=sign(load_conditions(imaxP,2))*pu;
        mux=max(abs(load_conditions(:,3)));
        muy=max(abs(load_conditions(:,4)));

        excentricity_x=abs(mux/pu);
        excentricity_y=abs(muy/pu);
        eccentricity_xy=[excentricity_x,excentricity_y];

        [Inertia_xy_modif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,...
            t1_var,eccentricity_xy,t2_var,pu,best_c,condition_cracking,E);
    
    end
end
if plotRebarDesign==1 && noptions>1
    diagramsFinalRebarCols(load_conditions,bestdiagrama,best_disposicion,...
                h,b,bestArrangement);
end