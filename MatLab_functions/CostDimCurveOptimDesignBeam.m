function [collectionDimBeams,collectionISRareaBeams,collectionISRbeams,...
          collectionEffBeams]=CostDimCurveOptimDesignBeam(span,wac,fc,...
          b_rec,h_rec,duct,pu_beams_steel,unit_cost_conc_beams,...
          load_conditions,fy)

%------------------------------------------------------------------------
% Syntax:
% [collectionDimBeams,collectionISRareaBeams,collectionISRbeams,...
%  collectionEffBeams]=CostDimCurveOptimDesignBeam(span,wac,fc,b_rec,...
%  h_rec,duct,pu_beams_steel,unit_cost_conc_beams,...
%  graphConvergencePlot,load_conditions,fy)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To compute a Cost-Dimension curve for reliable design options 
% of a beam structural element subject to pure-flxure
% 
% OUTPUT: collectionDimBeams:       Is the array containing the pair of
%                                   dimensions (b,h) for each potential
%                                   design option in the Cost-Curve
%
%         collectionISRareaBeams:   Is the array containing the ISR area
%                                   for each of the three cross-sections 
%                                   of each potential design option in the 
%                                   Cost-Dimension Curve
%
%         collectionISRbeams:       Is the array containing the ISR width 
%                                   values for each of the three cross-
%                                   sections of each potential design option 
%                                   in the Cost-Dimension Curve
%
%         collectionEffBeams:       Is the array containing the structural 
%                                   efficiencies for each of the three cross-
%                                   sections of each potential design option 
%                                   in the Cost-Dimension Curve
%
% INPUT:  span:                     is the span or length of the beam 
%                                   element (cm)
%
%         fc:                       is the concrete compression resistance 
%                                   to be used for the beam element f'c 
%                                   (Kg/cm2)
%
%         b_rec:                    Is the concrete cover along the width
%                                   dimension of the beam's cross-section 
%
%         h_rec:                    Is the concrete cover along the height
%                                   dimension of the beam's cross-section
%
%         duct:                     is the ductility demand level of design: 
%                                   1 - low ductility, 2 - medium ductility, 
%                                   3 - high ductility
%
%         pu_beams_steel:           is the unit cost of rebar assembly in
%                                   beams
%
%         fy:                       is the reinforcement steel yielding 
%                                   stress in (Kg/cm)
%
%         load_conditions:          are the load conditions for all three 
%                                   cross-sections: a vector consisting of
%                                   one row and four columns as: 
%                                   [1 Mu{left} Mu{mid} Mu{right}] in units 
%                                   Ton-m. The sign of each load should be
%                                   included
%
%         unit_cost_conc_beams:     is the unit cost of the concrete
%                                   (according to the f'c value) per unit
%                                   volume
%
%------------------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-14
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-04-15
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% Initializing dimensions
bmin=20; % from code (cm)
bi=bmin;
hi=bi;

CostISRVec=[];
CostConcreteVec=[];
TotalCostStrucVec=[];

% OPTIMAL DESIGN
factor_fc=0.85;
fdpc=factor_fc*fc;

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

collectionISRbeams=[];
collectionISRareaBeams=[];
collectionDimBeams=[];
collectionEffBeams=[];
bmax_condition=0;
while bmax_condition==0
    bp=bi-2*b_rec;

    %%%%%%%%%%%%%%%%%%%%%%%% ISR1t-Beams - SGD %%%%%%%%%%%%%%%%%%%%%%%%
    bh_condition=0;
    while bh_condition==0

        E=2.0e6; % Modulus of Elasticity of reinforcing steel
        %----------------------------------------------------------------%
        inertia_modif=zeros(1,3);

        Mr_3section=zeros(1,3);
        Ef_3section=zeros(1,3);
        collectionISRtension=[];
        
        %% ISR design for negative moment (right section)

        %%%%%%%%%%%%%%%%%%%%%%% Variable t range %%%%%%%%%%%%%%%%%%%%%%%%%

        d=hi-h_rec;
        tmin=(0.7*sqrt(fdpc)/fy*bi*d)/bp;
        if duct==1
            tmax=0.9*bi*d*fdpc/fy*(6000*betac)/((fy+6000))/bp;
        elseif duct==2
            tmax=(0.75*bi*d*fdpc/fy*(6000*betac))/((fy+6000))/bp;
        elseif duct==3
            tmax=(0.025*bi*d)/bp;
        end
        conditions_right=[load_conditions(1,1) load_conditions(1,4)];

        [cbest,bestMr,bestef03,best_Area03,tbest,hi]=SGD1tBeamsISR(bi,hi,duct,...
          b_rec,h_rec,fc,conditions_right,factor_fc,E,0);

        % Compute cost with the ISR
        bestCost=pu_beams_steel*best_Area03*wac*span;

        t2Best=[tbest,tmin];
        collectionISRtension=[collectionISRtension,tbest];

        Ef_3section(3)=bestef03;
        Mr_3section(3)=bestMr;
        h_right=hi;

        bestAreatDer_ten=t2Best(1)*bp;
        bestAreatDer_com=t2Best(2)*bp;

        areaISR_Der=[bestAreatDer_ten bestAreatDer_com];
        twidths_Der=[t2Best(1) t2Best(2)];
        
        unitCostDer=bestCost;

        bestAreatDer=bestAreatDer_com+bestAreatDer_ten;

        Inertia_modif=InertiaBeamCrackedSection(fc,E,areaISR_Der,bi,hi,h_rec);
        inertia_modif(3)=Inertia_modif;

        %% ISR design for positive moment (central section)

        d=hi-h_rec;
        tmin=(0.7*sqrt(fdpc)/fy*bi*d)/bp;

        conditions_central=[load_conditions(1,1) load_conditions(1,3)];

        [cbest,bestMr,bestef02,best_Area02,tbest,hi]=SGD1tBeamsISR(bi,hi,duct,...
            b_rec,h_rec,fc,conditions_central,factor_fc,E,0);

        % Compute cost with the ISR
        bestCost=pu_beams_steel*best_Area02*span*wac;

        t2Best=[tbest,tmin];
        collectionISRtension=[collectionISRtension,tbest];

        bestAreatCen_ten=t2Best(1)*bp;
        bestAreatCen_com=t2Best(2)*bp;

        twidths_Cen=[t2Best(1) t2Best(2)];
        areaISR_Cen=[bestAreatCen_ten bestAreatCen_com];

        Ef_3section(2)=bestef02;
        Mr_3section(2)=bestMr;
        h_central=hi;
        unitCostCentral=bestCost;

        bestAreatCentral=bestAreatCen_com+bestAreatCen_ten;

        %%% modified inertia (mid section)
        Inertia_modif=InertiaBeamCrackedSection(fc,E,areaISR_Cen,bi,hi,h_rec);
        inertia_modif(2)=Inertia_modif;

        %% Bar design for positive moment (left section)

        d=hi-h_rec;
        tmin=(0.7*sqrt(fdpc)/fy*bi*d)/bp;

        conditions_izq=[load_conditions(1,1) load_conditions(1,2)];

        [cbest,bestMr,bestef01,best_Area01,tbest,hi]=SGD1tBeamsISR(bi,hi,duct,...
                b_rec,h_rec,fc,conditions_izq,factor_fc,E,0);

        % Compute cost with the ISR
        bestCost=pu_beams_steel*best_Area01*span*wac;
        t2Best=[tbest,tmin];
        collectionISRtension=[collectionISRtension,tbest];

        bestAreatIzq_ten=t2Best(1)*bp;
        bestAreatIzq_com=t2Best(2)*bp;
        twidths_Izq=[t2Best(1) t2Best(2)];
        areaISR_Izq=[bestAreatIzq_ten bestAreatIzq_com];

        h_left=hi;

        Ef_3section(1)=bestef01;
        Mr_3section(1)=bestMr;
        unitCostIzq=bestCost;

        %%% modified inertia (left  section)
        Inertia_modif=InertiaBeamCrackedSection(fc,E,areaISR_Izq,bi,hi,h_rec);
        inertia_modif(1)=Inertia_modif;

        minAreaISR3sec=[bestAreatIzq_ten, bestAreatIzq_com, bestAreatCen_ten,...
                bestAreatCen_com bestAreatDer_ten bestAreatDer_com];

        bestCostISR=(unitCostIzq+unitCostCentral+unitCostDer)/3;
    
        % To make sure all three cross-section dimensions are the same (in case
        % they may suffer modifications through the optimization process)
        if h_left==h_right && h_right==h_central
            break;
        else
            h_sections=[h_left,h_central,h_right];
            hmax=max(h_sections);
            h_left=hmax;
            h_central=hmax;
            h_right=hmax;

            continue;
        end
    end
    %% Compute the concrete cost
    CostConcrete=unit_cost_conc_beams*span*bi*hi;

    %% Compute the total cost of construction: concrete + steel
    TotalCostStruc=CostConcrete+bestCostISR;
    
    % If at least in one of the three cross-sections the required steel
    % reinforcement area is equal to the maximum allowed, then the beam's
    % cross-section is regarded as being "Too small" and the dimensions to
    % be increased 
    if collectionISRtension(1)==tmax || ...
       collectionISRtension(2)==tmax || ...
       collectionISRtension(3)==tmax
   
        bi=bi+5;
        hi=bi;
        continue;
        
    % If at least in one of the three cross-sections the required steel
    % reinforcement area is equal to the minimum allowed, then the beam's
    % cross-section is regarded as being "Too big" and the dimensions to
    % be decreased 
    elseif collectionISRtension(1)==tmin || ...
           collectionISRtension(2)==tmin || ...
           collectionISRtension(3)==tmin
        
        break;
        
    % if Amin<As<Amax, then the option is reliable ...
    else
        TotalCostStrucVec=[TotalCostStrucVec;TotalCostStruc];
        CostConcreteVec=[CostConcreteVec;CostConcrete];
        CostISRVec=[CostISRVec;bestCostISR];
        
        collectionDimBeams=[collectionDimBeams;
                            bi hi];
                        
        collectionISRbeams=[collectionISRbeams;
                        twidths_Izq,twidths_Cen,twidths_Der];
                    
        collectionISRareaBeams=[collectionISRareaBeams;
                        areaISR_Izq,areaISR_Cen,areaISR_Der];
        
        collectionEffBeams=[collectionEffBeams;
                            Ef_3section];
        
        bi=bi+5;
        hi=bi;
        
        continue;
        
    end
end
%% Plotting Cost-Dimension Design Curve
figure(5)
plot(collectionDimBeams(:,1),CostISRVec,'b -o','LineWidth',1.5,...
    'MarkerFaceColor','blue')
hold on
xlabel('b (cm)')
ylabel('Construction Cost ($)')
title('Cost-Dimension curve of Optimal Beam Designs')
legend('Steel Assembly')

figure(5)
plot(collectionDimBeams(:,1),CostConcreteVec,'m -o','LineWidth',1.5,...
    'MarkerFaceColor','magenta','DisplayName','Concrete Cost')
hold on

figure(5)
plot(collectionDimBeams(:,1),TotalCostStrucVec,'k -o','LineWidth',1.5,...
    'MarkerFaceColor','black','DisplayName','Total Cost')
hold on