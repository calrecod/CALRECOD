function [Mr_col,h,Inertia_xy_modif,bestArea,lowestCost,ovMostEc,nvEc,...
    maxEfEc,bestArrangement,best_disposicion]=optimalrebarColsSym(b,h,...
    rec,act,E,npuntos,fdpc,beta1,pu_col_sym,load_conditions,condition_cracking,...
    plotRebarDesign)

%------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,Inertia_xy_modif,bestArea,lowestCost,ovMostEc,nvEc,...
%   maxEfEc,bestArrangement,best_disposicion]=optimalrebarColsSym(b,h,...
%   rec,act,E,npuntos,fdpc,beta1,pu_col_sym,load_conditions,condition_cracking,...
%   plotRebarDesign)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design composed of
%          only one type of rebar.
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
% INPUT:  rec:                  concrete cover of cross-section for both 
%                               axis direction: [coverX,coverY]
%
%         act:                  optima ISR reinforcement area
%
%         sepMin:               min separation of rebars constraint
%
%         E:                    Elasticity Modulus of reinforcement steel 
%                               E=2.0x10^6(Kg/cm2)$
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

pu=load_conditions(1,2);
mux=load_conditions(1,3);
muy=load_conditions(1,4);

Mr_col=zeros(1,2);

bp=b-2*rec(1);
hp=h-2*rec(2);

tableOptionsBar=zeros(7,8);
% Data base: type of rebar
tableOptionsBar(1:7,2)=[4;5;6;8;9;10;12];

% Database: unit weight kg/m of each type of rebar
tableOptionsBar(1:7,8)=[0.994;1.552;2.235;3.973;5.033;6.225;8.938];

% Database: UnitPrice PU $/kg for each type of rebar
% The following construction performances will be used (kg/jor)
% for each type of rebar:
% _____________________________________________
% #4 - 212, #5 - 216, #6 - 220
% #8 - 220, #9 - 220, #10 - 220, %12 - 220
%______________________________________________

nopciones=0;
while nopciones==0
    lowestCost=inf;
    maxEfEc=1;
    for i=1:7
        
        tableOptionsBar(i,1)=i;
        ov=tableOptionsBar(i,2);
        av=((ov/8)*2.54)^2*pi/4;
        dv=(ov/8*2.54);
        tableOptionsBar(i,3)=av;
        tableOptionsBar(i,4)=dv;

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

        sepMin=max([1.5*2.54, 3/2*dv]);
        tableOptionsBar(i,5)=nv;

        % There is a limit of the number of rebars that can be laid out
        % on each boundary of the cross-section, for each type of rebar
        maxVarillasSup=fix((bp+sepMin)/(sepMin+dv));
        maxVarillasCos=fix((hp+sepMin)/(sepMin+dv));

        tableOptionsBar(i,6)=maxVarillasSup;

        minVarillasSup=0.5*(nv-2*maxVarillasCos);
        if (minVarillasSup<2)
            minVarillasSup=2;
        end
        tableOptionsBar(i,7)=minVarillasSup;
        if (2*maxVarillasSup+2*maxVarillasCos<nv)
            tableOptionsBar(i,:)=0;
            continue;
        else
            
            costo=nv*av*0.0001*7800*pu_col_sym(i);
            
            for type=minVarillasSup:maxVarillasSup
                varSup=type;
                varCos=0.5*(nv-2*varSup);

                arregloVar=zeros(nv,1)+i;
                [disposicion_varillado]=RebarDisposition(b,...
                     h,rec,dv,nv,varCos,varSup);
                if nv~=length(disposicion_varillado)
                    break;
                end
                [diagrama,maxef,eficiencia,cxy]=diagramasDisposicion(ast,b,h,E,npuntos,...
                                        fdpc,nv,beta1,ov,av,disposicion_varillado,load_conditions);

                if costo<=lowestCost && maxef<maxEfEc
                    nvEc=nv;
                    lowestCost=costo;
                    ovMostEc=ov;
                    bestMrx=eficiencia(1,5);
                    bestMry=eficiencia(1,7);
                    dvEc=dv;
                    avEc=av;
                    maxEfEc=maxef;
                    vxEc=varSup;
                    vyEc=varCos;
                    best_disposicion=disposicion_varillado;
                    bestdiagrama=diagrama;
                    bestArrangement=arregloVar;
                    bestArea=ast;
                    best_c=cxy;
                end
            end
        end
    end

    Mr_col(1)=bestMrx;
    Mr_col(2)=bestMry;
    
    %%% MODIFIED INERTIA MOMENTUM THROUGH THE TRANSFORMED CRACKED
    %%% CROSS-SECTION .............................................

    %%%% Calculation of steel area on each cross-section's boundary
    area_vx=vxEc*avEc; 
    area_vy=vyEc*avEc;

    t1_var=area_vx/(b-2*rec(1));
    t2_var=area_vy/(h-2*rec(2));
    
    excentricity_x=abs(mux/pu)*100;
    excentricity_y=abs(muy/pu)*100;
    eccentricity_xy=[excentricity_x,excentricity_y];
    
    [Inertia_xy_modif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,t1_var,eccentricity_xy,...
                            t2_var,pu,best_c,condition_cracking,E);
        
    disponibilidadVar=[];
    nopciones=0;
    for i=1:7
        if (tableOptionsBar(i,1)~=0)
            nopciones=nopciones+1;
            disponibilidadVar=[disponibilidadVar;
                           tableOptionsBar(i,:)];
        end
    end
    
    if nopciones==0
        fprintf('\nDimensiones no aptas para columna por varillado\n');
        fprintf('Se aumentará el peralte\n');

        h=h+5;
    end
end
if plotRebarDesign==1
diagramsFinalRebarCols(load_conditions,bestdiagrama,best_disposicion,...
                h,b,bestArrangement);
end
