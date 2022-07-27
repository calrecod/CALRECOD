function [b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_x,...
        t_value_y,cxy]=isrColumns(pu_cols,height,b,h,rec,fy,fc,load_conditions,...
        ductility,optimaConvPlot,plotISRResults)

%------------------------------------------------------------------------
% Syntax:
% [b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_x,...
%  t_value_y,cxy]=isrColumns(pu_cols,height,b,h,rec,fy,fc,load_conditions,...
%  ductility,optimaConvPlot,plotISRResults) 
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal ISR for a column cross-section given 
% certain load conditions through the SGD method.
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
%         fc:               is the f'c value (Kg/cm2)
%
%         beta1:            is determined as established in ACI 318 code
%                           (see Documentation)
%
%         load_conditions:  are the load conditions applied to a cross-
%                           section: size = $[nloads,4]$ as: [load,Pu,Mx,My]
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

npuntos=30;
Ac_sec_elem=0;
iteracion=0;
while Ac_sec_elem==0
    
    iteracion=iteracion+1;
    %%%%%%%%%%%%%%-------------------------------------------%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Propiedades de materiales %%%%%%%%%%%%%%%%%%%%%

    E=2.0e6;
    tma=0.75;

    sepMin=3/2*tma*2.54;
    dimensionesColumna=[b h];
    
    bp=b-2*rec(1);
    hp=h-2*rec(2);

    betac=1.05-fc/1400;
    if (betac<0.65)
        betac=0.65;
    elseif (betac>0.85)
        betac=0.85;
    end
    fdpc=betac*fc;
    
    efmax=1.0;
    
    limEficienciaMenor=0.9;
    limEficienciaMayor=1.0;
    
    % initialize values for interation
    if ductility==1 || ductility==2
        tmin=0.01*b*h/(2*bp+2*hp);
        tmax=0.06*b*h/(2*bp+2*hp);
    elseif ductility==3
        tmin=0.01*b*h/(2*bp+2*hp);
        tmax=0.04*b*h/(2*bp+2*hp);
        
    end
    t=0.5;
    tk1=t;
    
    [Eftk1,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(tk1,dimensionesColumna,rec,fy,npuntos,...
            load_conditions,fdpc,E,betac);
    best_diagramaInteraccion=diagramaInteraccion;
    
    best_efk2_table=tablaEficiencias;
    dt=0.00001;
    tk1dt=tk1+dt;
    
    cxk1=cxy;
    
    [Eftk1dt,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(tk1dt,dimensionesColumna,rec,fy,npuntos,...
            load_conditions,fdpc,E,betac);

    
    %compute next first step to initialize loop
    dEfdt1=(Eftk1dt-Eftk1)/(tk1dt-tk1);
    
    % initial step-length 0.1
    alfa0=0.05;
    if Eftk1>limEficienciaMayor
        
        pk1=1; % search direction
    elseif Eftk1<limEficienciaMenor
        pk1=-1; % search direction
    else
        pk1=1;
    end
    tk2=tk1+alfa0*pk1;
    j=0;
    tk2_vector=[];
    Eftk2_vector=[];
    while (Eftk1>limEficienciaMayor || Eftk1<limEficienciaMenor)
        if tk2<=tmin
            fprintf('La sección es muy grande, se recomienda reducir las dimensiones\n');
            tk2=tmin;
            break;
        elseif tk2>=tmax
            fprintf('La sección es muy pequeña, se recomienda aumentar las dimensiones\n');
            tk2=tmax;
            break;
        end
        j=j+1;
        % Evaluate next function step
        [Eftk2,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(tk2,dimensionesColumna,rec,fy,npuntos,...
            load_conditions,fdpc,E,betac);
        
        cxk2=cxy;
        
        best_efk2_table=tablaEficiencias;
        best_diagramaInteraccion=diagramaInteraccion;
        
        tdtk2=tk2+dt;
        [Eftdtk2,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(tdtk2,dimensionesColumna,rec,fy,npuntos,...
            load_conditions,fdpc,E,betac);

        dEfdt2=-abs((Eftdtk2-Eftk2)/(tdtk2-tk2));
        
        %direction vector next step
        if Eftk1>limEficienciaMayor
            pk2=1;
            
            %step length
            alfak=abs(alfa0*(dEfdt2*pk2)/(dEfdt1*pk1));
        elseif Eftk1<limEficienciaMenor
            pk2=-1;
            
            %step length
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
        
        h=h+5;
        fprintf("\nLa sección de la columna es muy pequeño,\n se aumentará el peralte\n");
     
        continue;
    elseif tk2>tmin || tk2<tmax
        break;
    end
end
Ef_sec_col=Eftk1;
t_final_value=tk1;

ak1=2*bp*tk1+2*hp*tk1;

Ac_sec_elem=ak1;
cxy=cxk1;

%%%%%%%% Diagrama de interacción sección con perfil de acero %%%%%%%%%%%%%
if plotISRResults==1
    diagramISR(best_diagramaInteraccion,load_conditions);
end
t_value_x=t_final_value;
t_value_y=t_final_value*(h-2*rec(2))/(h-2*rec(2)-2*sepMin);

cost_elem_col=Ac_sec_elem*0.0001*height*0.01*7800*pu_cols;
Mr_col=[best_efk2_table(1,5) best_efk2_table(1,7)];

if optimaConvPlot==1
    figure(1)
    plot(tk2_vector,...
        Eftk2_vector,'r o','linewidth',0.1,'MarkerFaceColor','red' )
    hold on
    xlabel('t width (cm)')
    ylabel('Efficiency (%)')
    title('Search space for optimum width t')
end
