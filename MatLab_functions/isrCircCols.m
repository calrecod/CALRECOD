function [cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value,c]=isrCircCols...
    (pu_cols,height,wac,diam,rec,fy,fc,load_conditions,duct,optimaConvPlot,...
    plotISRResults)

%------------------------------------------------------------------------
% Syntax:
% [cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value,c]=isrCircCols...
% (pu_cols,height,wac,diam,rec,fy,fc,load_conditions,duct,optimaConvPlot,...
% plotISRResults)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal ISR through the Steepest Gradient 
% Descent method for a column of circular cross-section given certain 
% load conditions.
% 
% OUTPUT: cost_elem_col:    is the total construction cost of the element,
%                           considering both concrete and reinforcing steel 
%
%         Ac_sec_elem:      is the optimal reinforcement area
%
%         Ef_sec_col:       is the optimal structural efficiency for the
%                           cross-section
%
%         Mr_col:           is the resisting moment of the optimal 
%                           reinforced concrete column
%
%         t_value:          is the optimal ISR width t
%
%         c:                neutral axis depth of optimal design 
%                           corresponding to the critical load condition
%
% INPUT:  diam:              cross-section diameter 
%
%         fc:               is the f'c value
%
%         load_conditions:  are the load conditions applied to a cross-
%                           section: size = [nloads,3] as: [n-load,Pu,Mu]
%
%         duct:             parameter that indicates the ductility demand: 
%                           (1),(2),(3) for low, medium and high ductility
%
%         rec:              is the concrete cover for both axis directions:
%                           [coverX,coverY]
%
%         optimaConvPlot:   is the parameters that indicates if the optima
%                           convergence plot is required or not. Option are:
%                           (1) the plot is required, (2) they plot is not
%                           required
%
%         plotISRResults:   is the parameters that indicates if the ISR 
%                           interaction diagrams are required or not. 
%                           Options are: (1) they are required, (2) they 
%                           are not required
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

npdiag=30; % number of points to compute for the interaction diagram
E=fy/0.0021; % Modulus of Elasticity of the reinforcing steel

areaCol=diam^2*pi*0.25;

if fc<2000 % units: Kg,cm
    betac=1.05-fc/1400;
    if (betac<0.65)
        betac=0.65;
    elseif (betac>0.85)
        betac=0.85;
    end
else        % units: lb, in
    betac=0.85-0.05*(fc-4000)/1000;
    if betac<0.65
        betac=0.65;
    elseif betac>0.85
        betac=0.85;
    end
end
fdpc=0.85*fc;

limEffInf=0.9;
limEffSup=1.0;

% initialize values for iteration
if fc<2000 % units: Kg,cm
    if duct==1 || duct==2
        tmin=0.01*areaCol/(pi*diam);
        tmax=0.06*areaCol/(pi*diam);
    elseif duct==3
        tmin=0.01*areaCol/(pi*diam);
        tmax=0.04*areaCol/(pi*diam);
    end
    t=0.5; % initial ISR's width (cm)
else        % units: lb, in
    tmin=0.01*areaCol/(pi*diam);
    tmax=0.04*areaCol/(pi*diam);
    t=0.2; % initial ISR's width (in)
end

tk1=t;

[Eftk1,diagramInt,tabEficiencies,c]=diagCircColsISR(tk1,diam,rec,fy,npdiag,...
        load_conditions,fdpc,E,betac);
    
bestdiagramInt=diagramInt;

best_efk2_table=tabEficiencies;
dt=0.00001;
tk1dt=tk1+dt;

cxk1=c;

[Eftk1dt,diagramInt,tabEficiencies,c]=diagCircColsISR(tk1dt,diam,rec,fy,...
    npdiag,load_conditions,fdpc,E,betac);


% compute next first step to initialize loop
dEfdt1=(Eftk1dt-Eftk1)/(tk1dt-tk1);

% initial step-length
alfa0=0.05;
if Eftk1>limEffSup

    pk1=1; % search direction
elseif Eftk1<limEffInf
    pk1=-1; % search direction
else
    pk1=1;
end
tk2=tk1+alfa0*pk1;
j=0;
tk2_vector=[];
Eftk2_vector=[];
while (Eftk1>limEffSup || Eftk1<limEffInf)
 
    if tk2<=tmin
        fprintf('The cross-section is too big. \n');
        tk2=tmin;
        if Eftk1<limEffInf
            break;
        else
            continue;
        end
    elseif tk2>=tmax
        tk2=tmax;
        if Eftk1>limEffSup
            break;
        else
            continue;
        end
    end
    j=j+1;
    % Evaluate next function step
    [Eftk2,diagramInt,tabEficiencies,c]=diagCircColsISR(tk2,diam,rec,fy,...
        npdiag,load_conditions,fdpc,E,betac);

    cxk2=c;

    best_efk2_table=tabEficiencies;
    bestdiagramInt=diagramInt;

    tdtk2=tk2+dt;
    [Eftdtk2,diagramInt,tabEficiencies,c]=diagCircColsISR(tdtk2,diam,rec,...
        fy,npdiag,load_conditions,fdpc,E,betac);

    dEfdt2=-abs((Eftdtk2-Eftk2)/(tdtk2-tk2));

    % direction vector next step
    if Eftk1>limEffSup
        pk2=1;

        % step length
        alfak=abs(alfa0*(dEfdt2*pk2)/(dEfdt1*pk1));
    elseif Eftk1<limEffInf
        pk2=-1;

        % step length
        alfak=abs(alfa0*(dEfdt1*pk1)/(dEfdt2*pk2));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tk2_vector=[tk2_vector,tk2];
    Eftk2_vector=[Eftk2_vector,Eftk2];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cxk1=cxk2;
    tk1=tk2;
    tk2=tk1+alfak*pk2;

    pk1=pk2;

    % update values
    dEfdt1=dEfdt2;
    alfa0=alfak;

    Eftk1=Eftk2;
end

if tk2==tmax
    fprintf("\nThe cross-section diameter is too small!\n");
elseif tk2>tmin || tk2<tmax
    fprintf("\nThe reinforcement is balanced\n");
end

Ef_sec_col=Eftk1;
t_final_value=tk1;

ak1=(pi*(diam-2*rec))*tk1;

Ac_sec_elem=ak1;
c=cxk1;

%%%%%%%% Interaction Diagram with the optimum ISR %%%%%%%%%%%%%

if plotISRResults==1
    plotdiagramCircISR(bestdiagramInt,load_conditions)
end
t_value=t_final_value;

cost_elem_col=Ac_sec_elem*height*wac*pu_cols;
Mr_col=best_efk2_table(1,5);

%%%%%%%%%%%%%%%% Optimization convergence %%%%%%%%%%%%%%%%%%%%%%%%%
if optimaConvPlot==1
    figure(1)
    plot(tk2_vector,...
        Eftk2_vector,'r o','linewidth',0.1,'MarkerFaceColor','red' )
    hold on
    xlabel('t width')
    ylabel('Efficiency (%)')
    title('Optimum ISR width t of a circular columns cross-section')
end