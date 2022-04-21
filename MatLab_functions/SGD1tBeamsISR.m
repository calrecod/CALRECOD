function [c_best,bestMr,bestEf,best_area,tbest,h]=SGD1tBeamsISR(b,h,duct,...
        b_rec,h_rec,fc,load_conditions,factor_fc,E,graphConvergencePlot)

%------------------------------------------------------------------------
% Syntax:
% [c_best,bestMr,bestEf,best_area,tbest,h]=SGD1tBeamsISR(b,h,duct,...
%    b_rec,h_rec,fc,load_conditions,factor_fc,E,graphConvergencePlot)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement area for a given beam 
% cross-section with specified initially dimensions (b,h) through the SGD 
% method.
% 
% OUTPUT: c_best,bestMr,bestEf: The neutral axis depth for the optimal 
%                               design, the resistant bending moment for 
%                               the optimal design,
%
%         best_area,tbest:      The optimal reinforcement area, the optimal 
%                               t width of the ISR
%
%         h:                    The final cross-section height in case it 
%                               is modified from the given initial proposal 
%                               value
%
% INPUT:  load_conditions:      vector as [nload,Mu] size: nloads x 2
%
%         factor_fc:            is determined by de applicable design code. 
%                               The ACI 318-19 specifies it as 0.85
%
%         duct:                 is the ductility demand parameter, with 
%                               possible values of 1,2 or 3, for low 
%                               ductility, medium ductility or high 
%                               ductility respectively
%
%         h_rec,b_rec:          is the concrete cover along the height dimension
%                               and the width cross-section dimension, respectively
%                               (cm)
%
%         h,b:                  cross-section dimensions (cm)
%
%         E:                    is the Elasticity Modulus of reinforcing steel
%                               (Kg/cm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fy=4200;

if fc<=280 
    betac=0.85;
elseif fc>280
    betac=1.05-fc/1400;
    if betac<0.65
        betac=0.65;
    elseif betac>0.85
        betac=0.85;
    end
end

fdpc=fc*factor_fc;
                       
bp=b-2*b_rec;

d=h-h_rec;
if duct==1
    tmax=0.9*b*d*fdpc/fy*(6000*betac)/((fy+6000))/bp;
    tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
elseif duct==2
    tmax=(0.75*b*d*fdpc/fy*(6000*betac))/((fy+6000))/bp;
    tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
elseif duct==3
    tmax=(0.025*b*d)/bp;
    tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
end

% initialize values for interation
t=0.5;
t_comp=tmin;

t1=t_comp;
asc=t1*bp;
    
tmax_condition=0;
while tmax_condition==0
    tk1=t;
    t2=tk1;
    ast=t2*bp;
    [Eftk1,mrtk1,ck1]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,betac);
                                 
    dt=0.00001;
    tk1dt=tk1+dt;
    t2=tk1dt;
    ast=t2*bp;
    [Eftk1dt,mrtk1dt,c]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,betac);
                                 
    %compute next first step to initialize loop
    dEfdt1=(Eftk1dt-Eftk1)/(tk1dt-tk1);
    
    % initial step-length 1.0
    alfa0=0.05;
    if Eftk1>1.0
        
        pk1=1; % search direction
    elseif Eftk1<0.9
        pk1=-1; % search direction
    else
        pk1=1;
    end
    tk2=tk1+alfa0*pk1;
    
    tk2_vector=[];
    Eftk2_vector=[];
    j=0;
    while (Eftk1>1.0 || Eftk1<0.9)
        j=j+1;
        % Evaluate next function step
        
        t2=tk2;

        ast=t2*bp;
        [Eftk2,mrtk2,ck2]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,betac);
                                 
        
        tdtk2=tk2+dt;
        t2=tdtk2;
        ast=t2*bp;
        [Eftdtk2,mrdtk2,c]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,betac);
                                 
        dEfdt2=-abs((Eftdtk2-Eftk2)/(tdtk2-tk2));
        
        %direction vector next step
        if Eftk1>1.0
            pk2=1;
            
            %step length
            alfak=alfa0*(dEfdt2*pk2)/(dEfdt1*pk1);
        elseif Eftk1<0.9
            pk2=-1;
            
            %step length
            alfak=alfa0*(dEfdt1*pk1)/(dEfdt2*pk2);
        end
       
        tk2_vector=[tk2_vector,tk2];
        Eftk2_vector=[Eftk2_vector,Eftk2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ck1=ck2;
        tk1=tk2;
        
        tk2=abs(tk1+alfak*pk2);
        pk1=pk2;

        % update values
        dEfdt1=dEfdt2;
        alfa0=alfak;
        mrtk1=mrtk2;
        
        Eftk1=Eftk2;
        
        %termiinal condition
        if tk1>tmax
            tk1=tmax;
            break;
        elseif  tk1<tmin
            tk1=tmin;
            break;
        end
    end
    
    if tk1==tmax
        h=h+5;
        d=h-h_rec;
        
        % Update minimum and maximum area by code
        if duct==1
            tmax=(0.9*b*d*fdpc/fy*(6000*betac))/(fy+6000)/bp;
            tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
        elseif duct==2
            tmax=(0.75*b*d*fdpc/fy*(6000*betac))/(fy+6000)/bp;
            tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
        elseif duct==3
            tmax=(0.025*b*d)/bp;
            tmin=(0.7*sqrt(fdpc)/fy*b*d)/bp;
        end
    else
        break;
    end
end
ak1=bp*(tk1);
bestMr=mrtk1;
tbest=tk1;
c_best=ck1;

best_area=ak1+asc;
bestEf=Eftk1;

if graphConvergencePlot==1
    % Optima Convergence Graph____________________________________________
    figure(4)
    plot(tk2_vector,...
       Eftk2_vector,'b o','MarkerFaceColor','blue')
    hold on
    xlabel('t width (cm)')
    ylabel('Efficiency (%)')
    title('Search space for optimum width t')
    legend('Beam Cross-section Optimization')
end
    
    
        
