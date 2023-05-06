function [c_best,bestMr,bestEf,best_area,tbest]=SGD1tTBeamsISR(bp,ht,ba,...
    ha,Lb,duct,cover,fc,load_conditions,factor_fc,E,graphConvergencePlot)

%------------------------------------------------------------------------
% Syntax:
% [c_best,bestMr,bestEf,best_area,tbest]=SGD1tTBeamsISR(bp,ht,ba,...
%   ha,Lb,duct,cover,fc,load_conditions,factor_fc,E,graphConvergencePlot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement area for a given beam 
% of T cross-section with specified initially dimensions (b,h) through the 
% Steepest Gradient Descent method.
% 
% OUTPUT: c_best,bestMr,bestEf: The neutral axis depth for the optimal 
%                               design, the resistant bending moment for 
%                               the optimal design and the final structural
%                               efficiency for the optimal design
%
%         best_area,tbest:      The optimal reinforcement area, the optimal 
%                               t width of the ISR
%
% INPUT:  load_conditions:      vector as [nload,Mu] size: nloads x 2
%
%         factor_fc:            is determined by the applicable design 
%                               code. The ACI 318-19 specifies it as 0.85
%
%         duct:                 is the ductility demand parameter, with 
%                               possible values of 1,2 or 3, for low 
%                               ductility, medium ductility or high 
%                               ductility respectively
%
%         ba:                   is the effective flange width of the T-beam 
%                               cross-section
%
%         ht:                   is total height of the T-beam cross-section
%
%         bp:                   is the web width of the T-beam
%                               cross-section
%
%         ha:                   is the flange thickness of the T-beam
%                               cross-section
%
%         Lb:                   is the length of the beam element
%
%         E:                    is the Elasticity Modulus of reinforcing 
%                               steel
%
%         cover:                is the concrete cover for the reinforcement
%
%         graphConvergencePlot: is the parameter that indicates weather or
%                               not it is required to plot the 
%                               optimization convergence with the Steepest 
%                               Gradient Descent method 
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-20
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
fdpc=fc*factor_fc;
                       
bpp=bp-2*cover;

d=ht-cover;
fy=E*0.0021; % yield stress of the reinforcing steel

if fc<2000 % units Kg, cm
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
    % Max and min reinforcement area
    if duct==1
        tmax=0.9*(bp*d)*fdpc/fy*(6000*betac)/((fy+6000))/bpp;
        tmin=(0.7*sqrt(fdpc)/fy*(bp*d))/bpp;
    elseif duct==2
        tmax=(0.75*(bp*d)*fdpc/fy*(6000*betac))/((fy+6000))/bpp;
        tmin=(0.7*sqrt(fdpc)/fy*(bp*d))/bpp;
    elseif duct==3
        tmax=(0.025*(bp*d))/bpp;
        tmin=(0.7*sqrt(fdpc)/fy*(bp*d))/bpp;
    end
else % units lb, in
    betac=0.85-0.05*(fc-4000)/1000;
    if betac<0.65
        betac=0.65;
    elseif betac>0.85
        betac=0.85;
    end
    % Max and min reinforcement area
    tmax=200/fy*bp*d/bpp;
    tmin=3*sqrt(fc)/fy*bp*d/bpp;
    
end

% initialize values for interation
limInf=0.85; limSup=1.0; % max and min values of acceptable structural 
                         % efficiency 
t=0.5;
t_comp=tmin;

t1=t_comp;
asc=t1*bpp;

% begin optimization process
tk1=t;
t2=tk1;
ast=t2*bpp;
[Eftk1,mrtk1,ck1]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
                    ha,Lb,ast,asc,cover,betac);

dt=0.00001;
tk1dt=tk1+dt;
t2=tk1dt;
ast=t2*bpp;
[Eftk1dt,mrtk1dt,c]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
                    ha,Lb,ast,asc,cover,betac);

% compute next first step to initialize loop
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
while (Eftk1>limSup || Eftk1<limInf)
    j=j+1;
    % Evaluate next function step
    t2=tk2;

    ast=t2*bpp;
    [Eftk2,mrtk2,ck2]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
                    ha,Lb,ast,asc,cover,betac);


    tdtk2=tk2+dt;
    t2=tdtk2;
    ast=t2*bpp;
    [Eftdtk2,mrdtk2,c]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
                    ha,Lb,ast,asc,cover,betac);

    dEfdt2=-abs((Eftdtk2-Eftk2)/(tdtk2-tk2));

    % direction vector next step
    if Eftk1>1.0
        pk2=1;

        % step length
        alfak=alfa0*(dEfdt2*pk2)/(dEfdt1*pk1);
    elseif Eftk1<0.9
        pk2=-1;

        % step length
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

    % termiinal condition
    if tk1>tmax && Eftk1>=limSup
        tk1=tmax;
        break;
    elseif  tk1<tmin && Eftk1<=limInf
        tk1=tmin;
        break;
    end
end

if tk1==tmax
    disp('The cross-section dimensions are too small!')
end
ak1=bpp*(tk1);
bestMr=mrtk1;
tbest=tk1;
c_best=ck1;

best_area=ak1+asc;
bestEf=Eftk1;

% Optima Convergence Graph 
if graphConvergencePlot==1
    figure(4)
    plot(tk2_vector,...
       Eftk2_vector,'b o','MarkerFaceColor','blue')
    hold on
    xlabel('t width')
    ylabel('Efficiency (%)')
    title('Search space for optimum width t of a T beam cross-section')
end
