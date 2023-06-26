function [b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_x,...
        t_value_y,cxy]=isrColumns(pu_cols,height,b,h,rec,fy,fc,...
        load_conditions,wac,ductility,optimaConvPlot,plotISRResults)

%------------------------------------------------------------------------
% Syntax:
% [b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_x,...
%  t_value_y,cxy]=isrColumns(pu_cols,height,b,h,rec,fy,fc,...
%  load_conditions,wac,ductility,optimaConvPlot,plotISRResults) 
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal ISR for a column cross-section given 
% certain load conditions through the SGD method.
% 
% Note: The structural efficiency is determined with the Inverse Load
% method (Bresler's formula) and the Contour Load method. Thus, only one
% interaction diagram is computed for the whole set of given load
% combinations.
%
% OUTPUT: b,h:              are the final cross-section dimensions in case 
%                           of a need of modification to comply with the 
%                           restrictions criteria
%
%         cost_elem_col:    is the total construction cost of the element,
%                           considering both concrete and reinforcing steel 
%
%         Ac_sec_elem:      is the optimal reinforcement area
%
%         Ef_sec_col:       is the optimal structural efficiency for the
%                           cross-section
%
%         Mr_col:           are the resisting moments for both axis 
%                           directions of the column cross-section: [Mrx,Mry]
%
%         t_value_x:        is the optimal ISR width $t$ in the x-axis of 
%                           the cross-section
%
%         t_value_y:        is the optimal ISR width t in the y-axis of the
%                           cross-section
%
%         cxy:              neutral axis depth of optimal design for both 
%                           axis axis directions of the column's 
%                           cross-section, corresponding to the critical 
%                           load condition: [cx,cy]
%
% INPUT:  b,h:              cross-section dimensions of column (width 
%                             and height)
%
%         fc:               is the f'c value
%
%         beta1:            is determined as established in ACI 318 code
%                           (see Documentation)
%
%         load_conditions:  are the load conditions applied to a cross-
%                           section: size = [nloads,4] as: [load,Pu,Mx,My]
%
%         ductility:        parameter than indicates the ductility demand: 
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
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

load_conditions(:,2:4)=abs(load_conditions(:,2:4));
npuntos=30; % number of points to be computed for the interaction diagram
Ac_sec_elem=0;
iter=0;
while Ac_sec_elem==0
    %% ISR optimization process with the Steepest Gradient Descent method
    iter=iter+1;
    
    % Materials 
    Es=fy/0.0021; % Modulus of elasticity of reinforcing steel
    tma=0.75; % in
    if fc<2000 % units: (Kg,cm)
        sepMin=3/2*tma*2.54;
    else       % units: (lb,in)
        sepMin=3/2*tma;
    end
    dimensionsCol=[b h];
    
    bp=b-2*rec(1);
    hp=h-2*rec(2);
    
    fdpc=0.85*fc; % reduced concrete's compressive strength
    if fc<2000 % units: Kg,cm
        betac=1.05-fc/1400;
        if (betac<0.65)
            betac=0.65;
        elseif (betac>0.85)
            betac=0.85;
        end
    else       % units: lb,in
        betac=0.85-0.05*(fc-4000)/1000;
        if betac<0.65
            betac=0.65;
        elseif betac>0.85
            betac=0.85;
        end
    end
    % Range of acceptable structural efficiency
    limEffInf=0.9;
    limEffSup=1.0;
    
    % Max and min values of ISR's t width value for termination
    if fc<2000 % units: Kg,cm (Mexican NTC-2017)
        if ductility==1 || ductility==2
            tmin=0.01*b*h/(2*bp+2*hp);
            tmax=0.06*b*h/(2*bp+2*hp);
        elseif ductility==3
            tmin=0.01*b*h/(2*bp+2*hp);
            tmax=0.04*b*h/(2*bp+2*hp);
        end
    else       % units: lb,in (ACI 318)
        tmin=0.01*b*h/(2*bp+2*hp);
        tmax=0.04*b*h/(2*bp+2*hp);
    end
    % initialize values for interation
    if fc<2000 % units: kg,cm
        t=0.5;
    else       % units: lb,in
        t=0.2; 
    end
    tk1=t;
    
    [Eftk1,diagramInterac,tableEff,cxy]=widthEfficiencyCols(tk1,...
        dimensionsCol,rec,fy,npuntos,load_conditions,fdpc,Es,betac);
    
    best_diagramaInteraccion=diagramInterac;
    
    best_efk2_table=tableEff;
    dt=0.00001;
    tk1dt=tk1+dt;
    
    cxk1=cxy;
    
    [Eftk1dt,diagramInterac,tableEff,cxy]=widthEfficiencyCols(tk1dt,...
        dimensionsCol,rec,fy,npuntos,load_conditions,fdpc,Es,betac);

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
    subiter=0;
    tk1_vector=[tk1];
    Eftk1_vector=[Eftk1];
    
    while (Eftk1>limEffSup || Eftk1<limEffInf)
        subiter=subiter+1;
        % Evaluate next function step
        [Eftk2,diagramInterac,tableEff,cxy]=widthEfficiencyCols(tk2,...
            dimensionsCol,rec,fy,npuntos,load_conditions,fdpc,Es,betac);
        
        cxk2=cxy;
        
        best_efk2_table=tableEff;
        best_diagramaInteraccion=diagramInterac;
        
        tdtk2=tk2+dt;
        [Eftdtk2,diagramInterac,tableEff,cxy]=widthEfficiencyCols(tdtk2,...
          dimensionsCol,rec,fy,npuntos,load_conditions,fdpc,Es,betac);

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
        
        % update values
        cxk1=cxk2;
        tk1=tk2;
        
        tk2=tk1+alfak*pk2;
        
        pk1=pk2;

        dEfdt1=dEfdt2;
        alfa0=alfak;
        
        Eftk1=Eftk2;
        
        % Storing convergence values for plot
        tk1_vector=[tk1_vector,tk1];
        Eftk1_vector=[Eftk1_vector,Eftk1];
        
        % Termination conditions
        if tk1<tmin && Eftk1<limEffInf
            tk1=tmin;
            break;
        elseif tk1>tmax && Eftk1>limEffSup
            tk1=tmax;
            break;
        elseif subiter>100
            break;
        end
    end
    
    if tk1>=tmax
        fprintf("\nThe column cross-section is too small");
        if fc<2000 % units: (Kg,cm)
            h=h+5; % height is increased 5 cm
        else       % units: (lb,in)
            h=h+2; % height is increased 2 inches
        end
        continue;
    elseif tk1>=tmin && tk1<=tmax
        fprintf("\nThe column cross-section reinforcement is balanced");
        break;
    elseif tk1<=tmin
        fprintf("\nThe column cross-section is too big. The min");
        fprintf("\nreinforcement will be used");
        tk1=tmin;
        break;
    end
    
end
Ef_sec_col=Eftk1;
t_final_value=tk1;

ak1=2*bp*tk1+2*hp*tk1;

Ac_sec_elem=ak1;
cxy=cxk1;

%% Interaction diagrams with the optimal ISR
if plotISRResults==1
    diagramISR(best_diagramaInteraccion,load_conditions);
end
t_value_x=t_final_value;
t_value_y=t_final_value*(h-2*rec(2))/(h-2*rec(2)-2*sepMin);

cost_elem_col=Ac_sec_elem*wac*height*pu_cols;
Mr_col=[best_efk2_table(1,5) best_efk2_table(1,7)];

%% Optimization convergence resutls
if optimaConvPlot==1
    figure(1)
    plot(tk1_vector,...
        Eftk1_vector,'r o','linewidth',0.1,'MarkerFaceColor','red' )
    hold on
    xlabel('t width')
    ylabel('Efficiency (%)')
    title('Search space for optimum width t')
end
