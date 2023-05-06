function [Mr_col,Inertia,bestArea,bestCost,bestdiagram,bestnv,bestobar,...
    bestEf,bestType,bestDisposition,bestc]=optimalRebarCirc...
    (diam,rec,act,Es,npdiag,fdpc,pu_col_circ,load_conditions,...
    optionsBar,wac,height,plotRebarDesign)
               
%------------------------------------------------------------------------
% Syntax:
% [Mr_col,Inertia,bestArea,bestCost,bestdiagram,bestnv,bestobar,...
%  bestEf,bestType,bestDisposition,bestc]=optimalRebarCirc...
%  (diam,rec,act,Es,npdiag,fdpc,pu_col_circ,load_conditions,...
%  optionsBar,wac,height,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal symmetrical rebar design for a circular
% concrete column.
% 
% OUTPUT: bestEf:               is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         bestCost:             is the construction cost corresponding to
%                               the optimal rebar design
%
%         Mr_col:               is the resisting moment of the optimal 
%                               reinforced concrete column
%                               
%         bestArea:             is the optimal rebar area
%
%         bestdiagram:          is the interaction diagram corresponding to
%                               the optimal rebar design
%               
%         bestnv:               is the number of rebar of the optimal 
%                               rebar design
%
%         bestobar:             is the eighth-of-inch of the rebar 
%                               corresponding to the optimal design
%
%         bestType:             is the list of rebar types' index of each
%                               rebar used in the optimal design - size:
%                               [nbars,1] (a list of numbers from 1 to 7)
%
%         bestDisposition:      is an array containing the local 
%                               coordinates of position of each rebar over 
%                               the cross-section corresponding to the 
%                               optimal rebar design 
%
% INPUT:  diam:                 is the diameter of the circular column
%
%         rec:                  is the concrete cover
%
%         act:                  is the ISR area
%
%         Es:                   is the Modulus of Elasticity of the 
%                               reinforcing steel
%
%         load_conditions:      is the array containing the load conditions:
%                               in format [nload,Pu,Mu]
%
%         npdiag:               is the number of points for the 
%                               interaction diagrams
%
%         fdpc:                 is 0.85*f'c
%          
%         pu_col_circ:          is the rebar unit construction cost:
%                               vector size: [1,nCommercialRebars]
%
%         plotRebarDesign:      is the parameter that indicates if the 
%                               rebar design results are required to plot
%                               not. Option are: (1) do plot, (2) do not
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-03-19
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fc=fdpc/0.85;
if fc<2000 % units: Kg,cm
    beta1=1.05-fc/1400;
    if (beta1<0.65)
        beta1=0.65;
    elseif (beta1>0.85)
        beta1=0.85;
    end
else       % units: lb, in
    beta1=0.85-0.05*(fc-4000)/1000;
    if beta1<0.65
        beta1=0.65;
    elseif beta1>0.85
        beta1=0.85;
    end
end
fy=Es*0.0021; % yield stress of reinforcing steel

bestobar=[];

ntrebar=length(optionsBar(:,1));
noptions=0;
bestArea=inf;
bestEf=1;
for i=1:ntrebar

    ov=optionsBar(i,1);
    dv=optionsBar(i,2);
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
    if fc<2000 % units: (Kg,cm)
        sepMin=max([1.5*2.54, 3/2*dv]);
    else       % units: (lb,in)
        sepMin=max([1.5, 3/2*dv]);
    end
    % There is a limit of the number of rebars that can be laid out
    % on each boundary of the cross-section, for each type of rebar
    maxVarillas=fix((diam-2*rec)/(sepMin+dv))+1;
    
    if (maxVarillas<6)
        maxVarillas=6;
    end
    if (maxVarillas)<nv
        continue;
    else
        noptions=noptions+1;
        costo=nv*av*height*wac*pu_col_circ(i);

        arregloVar=zeros(nv,1)+i;
        [dispositionRebar]=RebarDispositionCirc(diam,rec,dv,nv);
            
        [maxef,diagram,effTable,c]=diagCircColsRebar(ast,dispositionRebar,...
        diam,rec,fy,npdiag,load_conditions,fdpc,av,Es,beta1);
        
        critical=1;
        while maxef~=effTable(critical,5)
            critical=critical+1;
        end
        if ast<=bestArea && maxef<bestEf
            bestnv=nv;
            bestCost=costo;
            bestobar=ov;
            bestcritical=critical;
            bestMr=effTable(bestcritical,4);
            bestEf=maxef;
            bestDisposition=dispositionRebar;
            bestdiagram=diagram;
            bestType=arregloVar;
            bestArea=ast;
            bestc=c;
        end
    end
end

if noptions==0
    fprintf('\nThe diameter of the column are not fit for rebar.\n');

    Mr_col=[];          bestobar=[];            bestc=[];
    Inertia=[];         bestType=[];
    bestArea=[];        bestdiagram=[];
    bestCost=[];        bestnv=[];
    bestEf=[];          bestDisposition=[];
    
else
    if isempty(bestobar)==0
        Mr_col=bestMr;
        Inertia=1/4*pi*(0.5*diam)^4;
    else
        fprintf('\nThe dimensions of the column are not fit for rebar.\n');

        Mr_col=[];          bestobar=[];            bestc=[];
        Inertia=[];         bestType=[];
        bestArea=[];        bestdiagram=[];
        bestCost=[];        bestnv=[];
        bestEf=[];          bestDisposition=[];

    end
end

if plotRebarDesign==1 && isempty(bestDisposition)==0 
    
    plotdiagramCircRebar(load_conditions,bestdiagram,bestDisposition,...
                diam,bestType);
end