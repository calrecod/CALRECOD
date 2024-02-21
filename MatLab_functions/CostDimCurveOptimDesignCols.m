function [collectionDimCols,collectionISRareaCols,collectionISRcols,...
          collectionEffCols]=CostDimCurveOptimDesignCols(height,wac,fc,rec,...
          ductility,pu_cols,unit_cost_conc_cols,load_conditions,fy)

%------------------------------------------------------------------------
% Syntax:
% [collectionDimCols,collectionISRareaCols,collectionISRcols,...
%  collectionEffCols]=CostDimCurveOptimDesignCols(height,wac,fc,rec,...
%  ductility,pu_cols,unit_cost_conc_cols,load_conditions,fy)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To compute a Cost-Dimension Curve for reliable design options 
% of a concrete structural element subject to biaxial compression loads
% 
% OUTPUT: collectionDimCols:        Is the array containing the pair of
%                                   dimensions (b,h) for each potential
%                                   design option in the Cost-Dimension 
%                                   Curve
%
%         collectionISRareaCols:    Is the array containing the ISR area
%                                   for each potential design option in the 
%                                   Cost-Dimension Curve
%
%         collectionISRcols:        Is the array containing the ISR width 
%                                   values for each potential design option 
%                                   in the Cost-Dimension Curve
%
%         collectionEffCols:        Is the array containing the structural 
%                                   efficiency of each potential design  
%                                   option in the Cost-Curve
%
% INPUT:  height:                   is the column's height dimension
%
%         fc:                       is the f'c value
%
%         load_conditions:          are the load conditions applied to a 
%                                   cross-section: size = [nloads,4] as:
%                                   [load,Pu,Mx,My]
%
%         ductility:                parameter that indicates the ductility 
%                                   demand: 
%                                   (1),(2),(3) for low, medium and 
%                                   high ductility
%
%         rec:                      is the concrete cover for each
%                                   cross-section direction: 
%
%                                           [cover-b,cover-h]
%
%         fy:                       is the yield stress of the reinforcing
%                                   steel
%
%         unit_cost_conc_cols:      is the unit cost of the concrete per 
%                                   unit volume (according to the f'c
%                                   value)
%
%         pu_cols:                  is the unit construction cost of
%                                   reinforcing steel assembly
%       
%------------------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-15
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%
%------------------------------------------------------------------------

collectionDimCols=[];
collectionISRareaCols=[];
collectionISRcols=[];
collectionEffCols=[];
collectionTotalCostCols=[];

CostConcreteVec=[];
TotalCostStrucVec=[];
CostElemColSteelVec=[];

hmin=25; % cm
bmin=25; % cm

bi=bmin;
hi=hmin;

dimax=0;
ncols=0;
while dimax==0
    Ac_sec_elem=0;
    iteration=0;
    while Ac_sec_elem==0

        iteration=iteration+1;

        E=2.0e6; % Modulus of Elasticity of reinforcing steel
        tma=0.75; % agregate size (inches)
        sepMin=3/2*tma*2.54; % minimum separation (cm)
        npdiag=20; % number of points for the interaction diagrams

        dimensionesColumna=[bi hi];

        bp=bi-2*rec(1);
        hp=hi-2*rec(2);

        betac=1.05-fc/1400;
        if (betac<0.65)
            betac=0.65;
        elseif (betac>0.85)
            betac=0.85;
        end
        fdpc=0.85*fc; % reduced concrete compressive strength

        limEffInf=0.9; % sought range of structural efficiency
        limEffSup=1.0;

        % initialize values for interation
        if ductility==1 || ductility==2
            tmin=0.01*bi*hi/(2*bp+2*hp);
            tmax=0.06*bi*hi/(2*bp+2*hp);
        elseif ductility==3
            tmin=0.01*bi*hi/(2*bp+2*hp);
            tmax=0.04*bi*hi/(2*bp+2*hp);
        end
        t=0.5; % cm
        tk1=t;

        [Eftk1,diagramaInteraccion,tablaEficiencias,cxy]=...
            widthEfficiencyCols(tk1,dimensionesColumna,rec,fy,npdiag,...
                load_conditions,fdpc,E,betac);
            
        best_diagramaInteraccion=diagramaInteraccion;

        best_efk2_table=tablaEficiencias;
        dt=0.00001;
        tk1dt=tk1+dt;

        cxk1=cxy;

        [Eftk1dt,diagramaInteraccion,tablaEficiencias,cxy]=...
            widthEfficiencyCols(tk1dt,dimensionesColumna,rec,fy,npdiag,...
                load_conditions,fdpc,E,betac);

        % compute next first step to initialize loop
        dEfdt1=(Eftk1dt-Eftk1)/(tk1dt-tk1);

        % initial step-length 0.05
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
        tk1_vector=[tk1];
        Eftk1_vector=[Eftk1];
        while (Eftk1>limEffSup || Eftk1<limEffInf)
            
            j=j+1;
            % Evaluate next function step
            [Eftk2,diagramaInteraccion,tablaEficiencias,cxy]=...
                widthEfficiencyCols(tk2,dimensionesColumna,rec,fy,npdiag,...
                load_conditions,fdpc,E,betac);

            cxk2=cxy;

            best_efk2_table=tablaEficiencias;
            best_diagramaInteraccion=diagramaInteraccion;

            tdtk2=tk2+dt;
            [Eftdtk2,diagramaInteraccion,tablaEficiencias,cxy]=...
                widthEfficiencyCols(tdtk2,dimensionesColumna,rec,fy,npdiag,...
                load_conditions,fdpc,E,betac);

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

            cxk1=cxk2;
            tk1=tk2;
            tk2=tk1+alfak*pk2;

            pk1=pk2;

            % update values
            dEfdt1=dEfdt2;
            alfa0=alfak;

            Eftk1=Eftk2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Storing convergence values for plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tk1_vector=[tk1_vector,tk2];
            Eftk1_vector=[Eftk1_vector,Eftk2];
            
            % Termination conditions
            if tk1<tmin && Eftk1<limEffInf

                tk1=tmin;
                break;
            elseif tk1>tmax && Eftk1>limEffSup

                tk1=tmax;
                break;
            end
        end

        if tk1>=tmax
            fprintf("\nThe column cross-section is too small");
            hi=hi+5;
            
            continue;
        elseif tk1<=tmin
            
            fprintf("\nThe column cross-section is too big. The min");
            fprintf("\nreinforcement will be used");
            tk1=tmin;
            break;
        elseif tk1>=tmin && tk1<=tmax
            fprintf("\nThe column cross-section reinforcement is balanced");
            break;
        end
        
    end
    Ef_sec_col=Eftk1;
    t_final_value=tk1;
    ISRwidth=tk1;
    
    ak1=2*bp*tk1+2*hp*tk1;

    Ac_sec_elem=ak1;
    cxy=cxk1;

    t_value_x=t_final_value;
    t_value_y=t_final_value*(hi-2*rec(2))/(hi-2*rec(2)-2*sepMin);

    %% Compute the steel assembly construction cost
    cost_elem_col_steel=Ac_sec_elem*height*wac*pu_cols;
    
    %% Compute the concrete cost
    CostConcrete=unit_cost_conc_cols*height*bi*hi;
    
    %% Total construction cost: concrete + reinforcement
    TotalCostStruc=CostConcrete+cost_elem_col_steel;
    
    %% Cost-dimensions curve termination conditions
    % If the required steel
    % reinforcement area is equal to the maximum allowed, then the column's
    % cross-section is regarded as being "Too small" and the dimensions to
    % be increased 
    
    if ISRwidth>=tmax
        bi=bi+5;
        hi=hmin;
        continue;
    % If the required steel
    % reinforcement area is equal to the minimum allowed, then the beam's
    % cross-section is regarded as being "Too big" and the dimensions to
    % be decreased 
    elseif ISRwidth<=tmin
        break;

    % if Amin<As<Amax, then the option is reliable ...
    else
        
        if ncols>=2 && hi==hmin
            if hi==collectionDimCols(ncols,2)
                break;
            end
        end
        ncols=ncols+1;
        CostElemColSteelVec=[CostElemColSteelVec; cost_elem_col_steel];
        CostConcreteVec=[CostConcreteVec; CostConcrete];
        TotalCostStrucVec=[TotalCostStrucVec; TotalCostStruc];
        
        collectionDimCols=[collectionDimCols;
                            bi,hi];

        collectionISRareaCols=[collectionISRareaCols;
                                ak1];
        collectionISRcols=[collectionISRcols;
                            tk1];
        collectionEffCols=[collectionEffCols;
                            Ef_sec_col];

        collectionTotalCostCols=[collectionTotalCostCols;
                                TotalCostStruc];
        bi=bi+5;
        hi=hmin;

        continue;
    end
end
%% Plotting Cost-Dimension Design Curve
figure(5)
plot(collectionDimCols(:,1),...
   CostElemColSteelVec,'b -o','LineWidth',1.5,'MarkerFaceColor','blue')
hold on
xlabel('b (cm)')
ylabel('Construction Cost ($)')
title('Cost-Dimension Curve of Optimal Column Designs')
legend('Reinforcing Steel')

figure(5)
plot(collectionDimCols(:,1),...
   CostConcreteVec,'m -o','LineWidth',1.5,'MarkerFaceColor','magenta',...
   'DisplayName','Concrete Cost')
hold on

figure(5)
plot(collectionDimCols(:,1),...
   TotalCostStrucVec,'k -o','LineWidth',1.5,'MarkerFaceColor','black',...
   'DisplayName','Total Cost')
hold on