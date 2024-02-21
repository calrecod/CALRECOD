function [sepbarsRestricSections,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
    barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
    barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
    barArrangementIzqComp,minAreaVar_3sec,Ef_elem_sec_t,bestCostSteel,ef_var,...
    minAreaVar_prom,Mr_3section]=beamsISR(pu_beams,span,wac,b,h,h_rec,...
    rebarAvailable,fc,fy,load_conditions,cols_sym_asym_isr,duct,b_rec,...
    rebarDesignPlots,graphConvergencePlot)

%------------------------------------------------------------------------
% Syntax:
% [sepbarsRestricSections,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
% barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
% barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
% barArrangementIzqComp,minAreaVar_3sec,Ef_elem_sec_t,bestCostVar,ef_var,...
% minAreaVar_prom,Mr_3section]=beamsISR(pu_beams,span,wac,b,h,h_rec_sections,...
% rebarAvailable,fc,fy,load_conditions,cols_sym_asym_isr,duct,b_rec,...
% rebarDesignPlots,graphConvergencePlot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: Tto design optimally all three cross-sections (left, middle and 
% right) along a beam element according to its moment distribution diagram,
% using the ISR analogy.
% 
% OUTPUT: epbarsRestric:        is the parameter that indicates if the rebar
%                               separation constraint of rebars in tension 
%                               for each of the beam's cross-section is 
%                               being complied: (1) means the restriction 
%                               is not being complied in the design, (0) 
%                               means the restriction is being complied
%
%         inertia_modif:        is the modified cross-section inertia based 
%                               on a cracking mechanism for each of the 
%                               designed cross-sections
%
%         b,h:                  are the final cross-section dimensions in
%                               case they suffered modifications after the
%                               optimal design process
%
%         dispositionBar_Der:   are the local coordinates of rebar
%                               disposition over the optimal designed right 
%                               cross-section 
%
%         dispositionBar_Center:are the local coordinates of rebar 
%                               disposition over the optimal designed 
%                               central cross-section 
%
%         dispositionBar_Izq:   are the local coordinates of rebar 
%                               disposition over the optimal designed left
%                               cross-section 
%
%         arrangement_t1:       are the list of rebar type transformed from
%                               the ISR in tension: a vector consisting of 
%                               one column of length $nbars$ in tension
%
%         barArrangementDerComp,
%         barArrangementDerTens:are the list of rebar type transformed from
%                               the ISR in compression and tension, 
%                               respectively, for the optimally designed 
%                               right cross-section: a vector consisting of
%                               one column of length $nbars$ in compression
%                               and tension
%
%         barArrangementCentralComp,
%         barArrangementCentralTens: are the list of rebar type transformed 
%                                    from the ISR in compression and tension,
%                                    respectively, for the optimally designed 
%                                    central cross-section: a vector
%                                    consisting of one column of length 
%                                    nbars in compression and tension
%
%         barArrangementIzqComp,
%         barArrangementIzqTens:are the list of rebar type transformed from
%                               the ISR in compression and tension, 
%                               respectively, for the optimally designed 
%                               left cross-section: a vector consisting of
%                               one column of length $nbars$ in compression
%                               and tension
%
%         ef_var:               is the optimal final structural efficiency
%                               for each of the three optimal designed beam 
%                               cross-sections considering the optimal rebar
%
%         Mr_3section:          is the optimal final bending resistance for
%                               each of the three optimal designed beam 
%                               cross-sections considering either the optimal
%                               rebar or the optimal ISR, according to the
%                               user preferences
%
%         bestCostVar:          is the total final cost of reinforcement 
%                               considering the three cross-section 
%                               reinforcement as an average rebar area along 
%                               the total span length of the beam element
%
%         minAreaVar_prom:      is the average rebar area of all three 
%                               cross-section (sum of steel in tension and
%                               compression)
%
% INPUT:  fc:                   is the concrete f'c
%
%         fy:                   is the reinforcement steel yielding stress 
%
%         load_conditions:      are the load conditions for all three 
%                               cross-sections: a vector consisting of one
%                               row and four columns as: 
%                               [1 Mu{left} Mu{mid} Mu{right}]. The sign of
%                               each load should be included
%
%         pu_beams:             is the unitary cost of rebar assembly in
%                               beams considering an average of assembly
%                               performance (it is considered that various
%                               types of rebars are placed simultaneously
%                               in the beam element along its length)
%
%         span:                 is the span length of the beam element
%
%         b,h:                  are the initial cross-section dimensions
%                               of the beam element
%
%         h_rec_sections:       are the concrete height cover for each of
%                               the three cross-sections as a vector of six
%                               columns and one row as:
%                               [cover_left_up,cover_left_low,cover_mid_up,
%                               cover_mid_low,cover_right_up,cover_right_low]
%
%         cols_sym_asym_isr:    is the optimal design option: either "ISR"
%                               (when no optimal rebar design is required,
%                               but only the ISR) or "Standard" (when a rebar
%                               optimal design is required)
%
%         duct:                 is the ductility demand level of design: 
%                               1 - low ductility, 2 - medium ductility, 
%                               3 - high ductility
%
%         b_rec:                is the lateral concrete cover of the beam 
%                               element as: [b_cover] (as it is the same 
%                               value for all three cross-sections)
%
%         plots:                is the parameter that indicates if the 
%                               plotting of results are required. Options 
%                               are: (1) they are required, (2) they are 
%                               not required
%
%         graphConvergencePlot: is the parameter that indicates if the 
%                               plotting of optima convergence is required.
%                               Options are: (1) they are required, 
%                                            (2) they are not required
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

factor_fc=0.85;
fdpc=factor_fc*fc;
bp=b-2*b_rec;
sepbarsRestricSections=[];
%%%%%%%%%%%%%%%%%%%%%%%% ISR1t-Beams - SGD %%%%%%%%%%%%%%%%%%%%%%%%
bh_condition=0;
while bh_condition==0
    
    E=fy/0.0021; % Modulus of Elasticity of reinforcing steel
    %------------------------------------------------------------------%
    inertia_modif=zeros(1,3);
    
    Mr_3section=zeros(1,3);
    
    %_________________________________________________________________
    % Bar design for negative moment (right section)..................
    %_________________________________________________________________
    
    %%%%%%%%%%%%%%%%%%%%%%% Variable range %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d=h-h_rec;
    if fc<2000 % units: kg,cm
        tmin=(0.7*sqrt(fc)/fy*b*d)/bp;
    else       % units: lb,in
        tmin=3*sqrt(fc)/fy*b*d/bp;
    end
    conditions_right=[load_conditions(1,1) load_conditions(1,4)];
    
    [cbest,bestMr,bestef03,best_Area03,tbest,h]=SGD1tBeamsISR(b,h,duct,...
      b_rec,h_rec,fc,conditions_right,factor_fc,E,graphConvergencePlot);
        
    % Compute cost with the ISR
    bestCost=pu_beams*best_Area03*wac*span;
    
    t2Best=[tbest,tmin];
    
    bestAreatDer_ten=t2Best(1)*bp;
    bestAreatDer_com=t2Best(2)*bp;
    
    if cols_sym_asym_isr~="ISR"
        
        [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,arregloVar_t1,...
            arregloVar_t2,ef,bestMr,area_var_t]=ISR1tRebarBeamsOptimization...
            (E,b,h,fy,fc,b_rec,h_rec,conditions_right,t2Best,pu_beams,wac,...
            span,rebarAvailable);

        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        barTypes1Right=arregloVar_t1;
        barTypes2Right=arregloVar_t2;
    else
        sepbarsRestric=0; % by default
        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        arregloVar_t1=[];
        arregloVar_t2=[];
        bestBarDisposition=[];
        area_var_t=[bestAreatDer_ten bestAreatDer_com];
        ef=0;
    end
    
    Mr_3section(3)=bestMr;
    h_right=h;
    bestAreaVarDer_ten=area_var_t(1);
    bestAreaVarDer_com=area_var_t(2);
    
    bestAreaVarDer=bestAreaVarDer_com+bestAreaVarDer_ten;
    
    unitCostDer=bestCost;
    bestEf_Der_var=ef;
    dispositionBar_Der=bestBarDisposition;

    barArrangementDerTens=arregloVar_t1;
    barArrangementDerComp=arregloVar_t2;

    %%% modified inertia (right section)
    
    Inertia_modif=InertiaBeamCrackedSection(fc,E,area_var_t,b,h,h_rec);
    inertia_modif(3)=Inertia_modif;

    %__________________________________________________________________
    % Bar design for positive moment (central section)................
    %__________________________________________________________________

    d=h-h_rec;
    if fc<2000 % units: kg,cm
        tmin=(0.7*sqrt(fc)/fy*b*d)/bp;
    else       % units: lb,in
        tmin=3*sqrt(fc)/fy*b*d/bp;
    end
    conditions_central=[load_conditions(1,1) load_conditions(1,3)];

    [cbest,bestMr,bestef02,best_Area02,tbest,h]=SGD1tBeamsISR(b,h,duct,...
        b_rec,h_rec,fc,conditions_central,factor_fc,E,graphConvergencePlot);
        
    % Compute cost with the ISR
    bestCost=pu_beams*best_Area02*wac*span;
    
    t2Best=[tbest,tmin];
    
    bestAreatCen_ten=t2Best(1)*bp;
    bestAreatCen_com=t2Best(2)*bp;
    
    if cols_sym_asym_isr~="ISR"
        
        [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,arregloVar_t1,...
         arregloVar_t2,ef,bestMr,area_var_t]=ISR1tRebarBeamsOptimization...
         (E,b,h,fy,fc,b_rec,h_rec,conditions_central,t2Best,pu_beams,wac,...
         span,rebarAvailable);

        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        barTypes1Mid=arregloVar_t1;
        barTypes2Mid=arregloVar_t2;
    else
        sepbarsRestric=0;
        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        arregloVar_t1=[];
        arregloVar_t2=[];
        bestBarDisposition=[];
        area_var_t=[bestAreatCen_ten bestAreatCen_com];
        ef=0;
    end
    
    Mr_3section(2)=bestMr;
    h_central=h;
    unitCostCentral=bestCost;
    bestAreaVarCentral_ten=area_var_t(1);
    bestAreaVarCentral_com=area_var_t(2);
    
    bestAreaVarCentral=bestAreaVarCentral_com+bestAreaVarCentral_ten;
    
    bestEf_Centrlal_var=ef;
    dispositionBar_Center=bestBarDisposition;

    barArrangementCentralTens=arregloVar_t1;
    barArrangementCentralComp=arregloVar_t2;

    %%% modified inertia (mid section)
    Inertia_modif=InertiaBeamCrackedSection(fc,E,area_var_t,b,h,h_rec);
    inertia_modif(2)=Inertia_modif;
    
    %__________________________________________________________________
    %%%%% Bar design for positive moment (left section)................
    %__________________________________________________________________
    
    d=h-h_rec;
    if fc<2000 % units: kg,cm
        tmin=(0.7*sqrt(fc)/fy*b*d)/bp;
    else       % units: lb,in
        tmin=3*sqrt(fc)/fy*b*d/bp;
    end
    
    conditions_izq=[load_conditions(1,1) load_conditions(1,2)];

    [cbest,bestMr,bestef01,best_Area01,tbest,h]=SGD1tBeamsISR(b,h,duct,...
            b_rec,h_rec,fc,conditions_izq,factor_fc,E,graphConvergencePlot);
        
    % Compute cost with the ISR
    bestCost=pu_beams*best_Area01*wac*span;
    
    t2Best=[tbest,tmin];
    bestAreatIzq_ten=t2Best(1)*bp;
    bestAreatIzq_com=t2Best(2)*bp;
    
    if cols_sym_asym_isr~="ISR"
        
        [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,arregloVar_t1,...
            arregloVar_t2,ef,bestMr,area_var_t]=ISR1tRebarBeamsOptimization...
            (E,b,h,fy,fc,b_rec,h_rec,conditions_izq,t2Best,pu_beams,wac,...
            span,rebarAvailable);

        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        barTypes1Left=arregloVar_t1;
        barTypes2Left=arregloVar_t2;
    else
        sepbarsRestric=0;
        sepbarsRestricSections=[sepbarsRestricSections,sepbarsRestric];
        arregloVar_t1=[];
        arregloVar_t2=[];
        bestBarDisposition=[];
        area_var_t=[bestAreatIzq_ten bestAreatIzq_com];
        ef=0;
    end
    h_left=h;

    Mr_3section(1)=bestMr;
    unitCostIzq=bestCost;
    bestEf_Izq_var=ef;
    
    bestAreaVarIzq_ten=area_var_t(1);
    bestAreaVarIzq_com=area_var_t(2);
    
    bestAreaVarIzq=bestAreaVarIzq_com+bestAreaVarIzq_ten;
    
    dispositionBar_Izq=bestBarDisposition;

    barArrangementIzqTens=arregloVar_t1;
    barArrangementIzqComp=arregloVar_t2;

    %%% modified inertia (left  section)
    Inertia_modif=InertiaBeamCrackedSection(fc,E,area_var_t,b,h,h_rec);
    inertia_modif(1)=Inertia_modif;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    minAreaVar_prom=(bestAreaVarIzq+bestAreaVarCentral+bestAreaVarDer)/3;
    
    if cols_sym_asym_isr~="ISR"
        minAreaVar_3sec=[bestAreaVarIzq_ten,bestAreaVarIzq_com,bestAreaVarCentral_ten,...
            bestAreaVarCentral_com, bestAreaVarDer_ten, bestAreaVarDer_com];
    
    else
        minAreaVar_3sec=[bestAreatIzq_ten, bestAreatIzq_com, bestAreatCen_ten,...
            bestAreatCen_com bestAreatDer_ten bestAreatDer_com];
        
    end
    
    bestCostSteel=(unitCostIzq+unitCostCentral+unitCostDer)/3;
    
    ef_var=[bestEf_Izq_var, bestEf_Centrlal_var, bestEf_Der_var];
    Ef_elem_sec_t=[bestef01 bestef02 bestef03];
    
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

% Plotting reinforced cross-section _______________________________

if rebarDesignPlots==1
    MidLeftRightBeamReinforcedSections(h,b,dispositionBar_Center,...
        dispositionBar_Izq,dispositionBar_Der,barTypes1Mid,barTypes2Mid,...
        barTypes1Left,barTypes2Left,barTypes1Right,barTypes2Right);
end
