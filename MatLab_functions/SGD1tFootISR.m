function [pbest,bestEf,bestMr,best_area]=SGD1tFootISR(b,h,rec,fdpc,fy,...
            steelAreaRange,betac,axis,mu_real_axis,plotOptimConv)

%------------------------------------------------------------------------
% Syntax:
% [pbest,bestEf,bestMr,best_area]=SGD1tFootISR(b,h,rec,fdpc,fy,...
%            steelAreaRange,betac,axis,mu_real_axis,plotOptimConv)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine the optimal steel area reinforcement of a footing 
% transversal cross-section subject to uni-axial flexure
% 
% OUTPUT: pbest:                is the optimal average reinforcement area
%                               percentage
%
%         bestEf:               is the structural efficiency of the optimal
%                               reinforcement option
%
%         bestMr:               is the resisting bending moment of the 
%                               optimal reinforced cross-section
%
%         best_area:            is the optima reinforcing area of the 
%                               transversal cross-section 
%
% INPUT:  h:                    is the cross-section height
%
%         b:                    is the cross-section width dimension
%
%         fdpc:                 is the reduced f'c as fdcp=0.85f'c according
%                               to code
%
%         rec                   is the concrete cover
%
%         mu_real_axis:         is the effective flexure distributed to the
%                               transversal cross-section from the contact
%                               soil pressures
%
%         fy:                   is the yield stress of reinforcing steel
%
%         axis:                 is the footing axis direction of analysis:
%                               (1) represents the axis direction in which 
%                               the dimension L is the width of the 
%                               transversal cross-section, for (2) B is the
%                               width of the transversal cross-section, in
%                               its reference system (see Documentation)
%
%         betac:                is determined as specified in code (see 
%                               Documentation)
%
%         plotOptimConv:        is the parameter that indicates if the 
%                               optima ISR convergence plot is required or
%                               not. Options are: (1) they are required, 
%                               (2) they are not required
%
%         steelAreaRange:       are the min and max percentage reinforcing
%                               area quantities for the cross-section
%                               according to ACI code (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%-------- geometry ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxEff=0.99;
minEff=0.8;

bp=b-2*rec;

if axis==1
    d=h-rec;
elseif axis==2
    d=h-rec-0.7*rec;
end
tmax=steelAreaRange(2)*b*(d)/bp;
tmin=steelAreaRange(1)*b*(d)/bp;

% initialize values for interation
t=0.5;
tk1=t;

[Eftk1,mrtk1]=EvaluateISR1tFoot(tk1,b,h,fy,fdpc,rec,betac,axis,mu_real_axis);

dt=0.00001;
tk1dt=tk1+dt;
[Eftk1dt,mrtk1dt]=EvaluateISR1tFoot(tk1dt,b,h,fy,fdpc,rec,betac,axis,mu_real_axis);

%compute next first step to initialize loop
dEfdt1=(Eftk1dt-Eftk1)/(tk1dt-tk1);

% initial step-length 1.0
alfa0=0.05;

if Eftk1>1.0

    pk1=1; % search direction
elseif Eftk1<0.8
    pk1=-1; % search direction
else
    pk1=1;
end
tk2=tk1+alfa0*pk1;
j=0;
tk2_vector=[];
Eftk2_vector=[];
while (Eftk1>maxEff || Eftk1<minEff)
    j=j+1;
    % Evaluate next function step
    [Eftk2,mrtk2]=EvaluateISR1tFoot(tk2,b,h,fy,fdpc,rec,betac,axis,mu_real_axis);

    tk2_vector=[tk2_vector,tk2];
    Eftk2_vector=[Eftk2_vector,Eftk2];

    tdtk2=tk2+dt;
    [Eftdtk2,mrdtk2]=EvaluateISR1tFoot(tdtk2,b,h,fy,fdpc,rec,betac,axis,mu_real_axis);

    dEfdt2=-abs((Eftdtk2-Eftk2)/(tdtk2-tk2));

    %direction vector next step
    if Eftk1>maxEff
        pk2=1;

        %step length
        alfak=alfa0*(dEfdt2*pk2)/(dEfdt1*pk1);
    elseif Eftk1<minEff
        pk2=-1;

        %step length
        alfak=alfa0*(dEfdt1*pk1)/(dEfdt2*pk2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tk1=tk2;
    tk2=abs(tk1+alfak*pk2);
    pk1=pk2;

    % update values
    dEfdt1=dEfdt2;
    alfa0=alfak;
    mrtk1=mrtk2;

    Eftk1=Eftk2;

    % termiinal condition

    if tk1>tmax
        continue;
    elseif tk1<tmin
        tk1=tmin;
        break;
    end
end

if plotOptimConv==1
    if axis==1
        figure(5)
        plot(tk2_vector,...
           Eftk2_vector,'r o','linewidth',0.1,'MarkerFaceColor','red')
        hold on
        xlabel('t width')
        ylabel('Efficiency (%)')
        title('Search space for optimum width t')
        legend('Section L')

    elseif axis==2
        figure(5)
        plot(tk2_vector,...
           Eftk2_vector,'b o','linewidth',0.1,'MarkerFaceColor','blue',...
           'DisplayName','Section B')
        hold on
        xlabel('t width')
        ylabel('Efficiency (%)')
        title('Search space for optimum width t')
    end
end
ak1=bp*(tk1);
bestMr=mrtk1;
tbest=tk1;
best_area=ak1;
bestEf=Eftk1;

pbest=ak1/(b*d);