function [hmodif,mu_axis,barDispositionFootings,arrangement_bar_footings,...
    nbars_footings,AcBar,bestCost_elem,list_ef_footings,list_mr_footings]=...
    isrFootings(pu_steel_footings,h,be,le,rec,fc,fy,load_conditions,dimCol,...
    RebarAvailable,cols_sym_asym_isr,ductility,optimConv,PlotRebarDesign,...
    typeFooting,sepMin_01,wac)

%------------------------------------------------------------------------
% Syntax:
% [hmodif,mu_axis,barDispositionFootings,arrangement_bar_footings,...
%   nbars_footings,AcBar,bestCost_elem,list_ef_footings,list_mr_footings]=...
%   isrFootings(pu_steel_footings,h,be,le,rec,fc,fy,load_conditions,dimCol,...
%   RebarAvailable,cols_sym_asym_isr,ductility,optimConv,PlotRebarDesign,...
%   typeFooting,sepMin_01,wac)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal rebar option of both transversal cross
% sections of an isolated footing subject to eccentric biaxial loads.
% 
% OUTPUT: hmodif:                   is the final modified design of the 
%                                   footing height
%
%         mu_axis:                  are the effective distributed bending 
%                                   moments for both transversal cross
%                                   section of the footing in analysis
%
%         barDispositionFootings:   is the collection of local rebar 
%                                   coordinates of both transversal cross
%                                   section of the footing. Size=[nbars,2]
%
%         arrangement_bar_footings: is the list that contains the rebar 
%                                   types of all rebars of the optimal 
%                                   designed footing. A number between 1-7
%                                   by default, comprising the 7 rebar types 
%                                   commercially available by default. 
%
%         nbars_footings:           is the list containing the total number 
%                                   of rebars used both in tension and 
%                                   compression for both transversal cross
%                                   sections of the footing
%
%         AcBar:                    is the list containing the quantity of
%                                   reinforcement area used for both 
%                                   transversal cross-sections as the sum
%                                   of the area in compression (or 
%                                   temperature) and tension
%         
%         bestCost:                 is the total assembly cost of the 
%                                   optimal reinforcement option, either 
%                                   with an ISR or with rebars
%
%         list_ef_footings:         is the list of final structural 
%                                   efficiencies for both transversal cross
%                                   sections of the footing. Size = [1,2]
%
%         list_mr_footings:         is the list of final resistant bending
%                                   moments for both of the transversal cross
%                                   sections. Size = [1,2]
%
% INPUT:  be,le:                    are the initial given transversal cross
%                                   section width dimensions
%
%         h:                        is the initial given footing height 
%                                   dimension
%
%         pu_steel_footings:        is the unit construction assembly cost 
%                                   of rebars in isolated footings: units 
%                                   $/weight. The unit cost is considered using 
%                                   an average value of all rebar types's 
%                                   assembly performances (assuming that in
%                                   an isolated footing as much as four 
%                                   different types of rebar may be placed)
%
%         fc:                       is the minimum steel area reinforcement
%                                   by temperature over the zone in 
%                                   compression of the transversal cross
%                                   section (upper boundary)
%
%         fy:                       is the yield stress of the reinforcing
%                                   steel
%
%         load_conditions:          are the eccentric biaxial load conditions
%                                   applied to the footing through the column
%                                   that supports
%
%         rec:                      is the concrete cover
%
%         dimCol:                   are the columns cross-section dimensions
%                                   that the footing supports
%
%         cols_sym_asym_isr:        is the parameter that indicates if only
%                                   an ISR optimal design is required or an
%                                   optimal rebar optimization design process.
%                                   Options are ''ISR'' or Symmetric
%
%         ductility:                is the ductility parameter that indicates
%                                   the level of ductility required for the
%                                   transversal cross-sections of the footing
%                                   through the max-min quantity of steel 
%                                   area. (1) low ductility, (2) medium 
%                                   ductility, (3) high ductility
%
%         optimConv:                is the parameters that indicates if the
%                                   optima ISR convergence plots are required
%                                   or not. Options are: (1) they required,
%                                   (2) they are not
%
%         PlotRebarDesign:          is the parameters that indicates if the
%                                   optima rebar design plots are required
%                                   or not. Options are: (1) they required,
%                                   (2) they are not
%
%         typeFooting:              1 - Standard Isolated footing
%                                   2 - Bordering Isolated footing
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%% ISR1t-Beams - SGD %%%%%%%%%%%%%%%%%%%%%%%%%
            
if fc<=280 
    beta1=0.85;
elseif fc>280
    beta1=1.05-fc/1400;
    if beta1<0.65
        beta1=0.65;
    end
end

d=h-rec; % effective footing's height
%% Distribution of contact soil pressures
[qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_conditions,be,le,...
                                              typeFooting,dimCol);
%% Shear design and revision
pu=load_conditions(1,2);
% Punching shear and shear by flexure
[d,qmax]=shearFootings(be,le,qprom,dimCol,pu,d,fc,typeFooting);

hmodif=d+rec; % new footing's height dimension (in case it was 
              % required by the shear design)
                  
fdpc=0.85*fc;
niter=0;    
maxiter=10;
amax_condition=0;
while amax_condition==0 % This loop stops until the max/min rebar
                        % area restriction is complied for both footing's
                        % cross-sections
    niter=niter+1; % counting number of iterations
    barDispositionFootings=[];
    arrangement_bar_footings=[];
    nbars_footings=zeros(1,4);
    bestCost_elem=0;
    AcBar=[];

    dim_zap=[be le];
    
    %%% Note: the analysis starts with the local Y-axis
    mu_axis=zeros(1,2);

    %% Reinforcement along Be - (section Le)
    eje=1;
    d=hmodif-rec;

    list_ef_footings=zeros(1,2);
    Ac_elem=[];
    
    % ------------------------------------------------------------------
    % MOMENT DISTRIBUTION 
    % ------------------------------------------------------------------
    if typeFooting==1
        bp=(be-dimCol(2)-d)*0.5;
    elseif typeFooting==2 || typeFooting==3
        bp=(be-dimCol(2)-d);
    end
    [m]=MomentDistributionFootings(qmax,bp,le);
    
    mu_axis(1)=m;

    pmin=0.0021;
    
    % Max reinforcing area (by code)
    if ductility==1 || ductility==2
        pmax=(0.9*le*d*fdpc/fy*(6000*beta1)/(fy+6000))/(le*d);
    elseif ductility==3
        pmax=(0.75*le*d*fdpc/fy*(6000*beta1)/(fy+6000))/(le*d);
    end

    variablesRange=[pmin pmax];
    [pbest,bestEf,bestMr,best_area]=SGD1tFootISR(le,hmodif,...
                    rec,fdpc,fy,variablesRange,beta1,eje,m,optimConv);
    p=pbest;
    if p<pmin
        p=pmin;
        amax_condition=1;
    elseif p>=pmax
        p=pmax;
        amax_condition=0;
        d=d+5;
        hmodif=d+rec;
        continue; % pass to the next while loop
    elseif p<=pmax && p>=pmin
        amax_condition=1;
    end
    
    ac=p*le*d;

    if cols_sym_asym_isr~="ISR" % When it is required to design the rebar
        % Reinforcing steel in tension 
        [le,ac_tension,nv,s,arreglo_t1]=...
         RebarOptionsFootings(ac,le,RebarAvailable,sepMin_01);
        
        if ac_tension==0
            d=d+5;
            hmodif=d+rec;
            continue;
        end
        [maxef,mr]=EfcriticalFootings(le,fdpc,...
            ac_tension,d,fy,m);

        list_ef_footings(1)=maxef;    % section Le
        list_mr_footings(1)=mr;
    
        % Minimum reinforcing area in compression by temperature
        acmin=pmin*le*d;

        [le,ac_comp,nv_comp,s,arreglo_t2]=...
            RebarOptionsFootings(acmin,le,RebarAvailable,sepMin_01);

        nv_section=[nv,nv_comp];

        if be>le
            largerDim=2; %%% axis 2
        elseif le>be
            largerDim=1; %%% axis 1
        else
            largerDim=0; % square footings
        end

        areaBar=ac_tension+ac_comp;
        AcBar=[AcBar,areaBar];

        if largerDim==0
            [dispositionRebar]=dispositionRebarSquareFootings...
                (le,hmodif,rec,RebarAvailable,nv_section,arreglo_t1,...
                arreglo_t2,eje);
        else
            if eje==largerDim
                [dispositionRebar,arreglo_t1]=...
                 dispositionRebarRectangularFootings(le,hmodif,rec,...
                 RebarAvailable,nv_section,arreglo_t1,arreglo_t2,eje,...
                 largerDim,dim_zap);
            else
                [dispositionRebar]=dispositionRebarSquareFootings(le,...
                hmodif,rec,RebarAvailable,nv_section,arreglo_t1,arreglo_t2,eje);
            end
        end
        
        barDispositionFootings=[barDispositionFootings;
                                dispositionRebar];

        arrangement_bar_footings=[arrangement_bar_footings;
                                        arreglo_t1;
                                       arreglo_t2];   

        nbars_footings(1,2*eje-1)=length(arreglo_t1);
        nbars_footings(1,2*eje)=length(arreglo_t2);

        bestCost_axis = EvaluateCostRebarFoot(be,RebarAvailable,arreglo_t1,...
            arreglo_t2,pu_steel_footings,wac);
        
        bestCost_elem=bestCost_elem+bestCost_axis;
        
        dispositionRebarFinal2=dispositionRebar;
        arreglo_t1L=arreglo_t1;
        arreglo_t2L=arreglo_t2;
        
    else
        acmin=pmin*le*d;
        areaBar=ac+acmin;
        AcBar=[AcBar,areaBar];
        
        list_ef_footings(1)=bestEf;    % section Le
        list_mr_footings(1)=bestMr;
        
        barDispositionFootings=[];

        arrangement_bar_footings=[];
            
        arreglo_t1=[];
        arreglo_t2=[];
        
        nbars_footings(1,2*eje-1)=length(arreglo_t1);
        nbars_footings(1,2*eje)=length(arreglo_t2);
        
        bestCost_axis = EvaluateCostISRFoot(be,rec,ac,acmin,...
            pu_steel_footings,wac);
                    
        bestCost_elem=bestCost_elem+bestCost_axis;

    end
    
    %% Reinforcement along le (section Be)
    eje=2;
    d=hmodif-rec-12/8*2.54;
    
    % -------------------------------------------------------------------
    % MOMENT DISTRIBUTION 
    % -------------------------------------------------------------------
    
    if typeFooting==1 || typeFooting==2
        lp=(le-dimCol(1)-d)*0.5;
    elseif typeFooting==3
        lp=(le-dimCol(1)-d);
    end
    [m]=MomentDistributionFootings(qmax,lp,be);

    mu_axis(2)=m;
    pmin=0.0021;

    if ductility==1 || ductility==2
        pmax=(0.9*be*d*fdpc/fy*(6000*beta1)/(fy+6000))/(be*d);
    elseif ductility==3
        pmax=(0.75*be*d*fdpc/fy*(6000*beta1)/(fy+6000))/(be*d);
    end

    variablesRange=[pmin pmax];
    [pbest,bestEf,bestMr,best_area]=SGD1tFootISR(be,hmodif,...
                    rec,fdpc,fy,variablesRange,beta1,eje,m,optimConv);
    p=pbest;
    if p<pmin
        p=pmin;
        amax_condition=1;
    elseif p>pmax
        p=pmax;
        amax_condition=0;
        d=d+5+12/8*2.54;
        hmodif=d+rec;
        continue; % pass to the next while loop
    elseif p<=pmax && p>=pmin
        amax_condition=1;
    end
    
    ac=p*be*d;
    if cols_sym_asym_isr~="ISR"  % When the design of rebar is required
        [be,ac_tension,nv_tension,s,arreglo_t1]=...
         RebarOptionsFootings(ac,be,RebarAvailable,sepMin_01);
        if ac_tension==0
            d=d+5+12/8*2.54;
            hmodif=d+rec;
            continue;
        end
        [maxef,mr]=EfcriticalFootings(be,fdpc,...
            ac_tension,d,fy,m);

        list_ef_footings(2)=maxef;
        list_mr_footings(2)=mr;

        % Minimum reinforcing steel by temperature in compression
        acmin=pmin*be*d;
        [be,a_comp,nv_comp,s,arreglo_t2]=...
            RebarOptionsFootings(acmin,be,RebarAvailable,sepMin_01);

        if be>le
            largerDim=2; %%% axis 2
        elseif le>be
            largerDim=1; %%% axis 1
        else
            largerDim=0; % square footings
        end

        areaBar=ac_tension+a_comp;
        AcBar=[AcBar,areaBar];
        nv_section=[nv_tension,nv_comp];

        if largerDim==0
            [dispositionRebar]=dispositionRebarSquareFootings(be,hmodif,rec,...
                     RebarAvailable,nv_section,arreglo_t1,arreglo_t2,eje);
        else
            if largerDim==eje
                [dispositionRebar,arreglo_t1]=...
                dispositionRebarRectangularFootings(be,hmodif,rec,...
                RebarAvailable,nv_section,arreglo_t1,arreglo_t2,eje,...
                largerDim,dim_zap);
            else
                [dispositionRebar]=dispositionRebarSquareFootings(be,...
                hmodif,rec,RebarAvailable,nv_section,arreglo_t1,...
                arreglo_t2,eje);
            end
        end
        
        barDispositionFootings=[barDispositionFootings;
                                dispositionRebar];

        arrangement_bar_footings=[arrangement_bar_footings;
                                        arreglo_t1;
                                       arreglo_t2];

        nbars_footings(1,2*eje-1)=length(arreglo_t1);
        nbars_footings(1,2*eje)=length(arreglo_t2);

        bestCost_axis = EvaluateCostRebarFoot(le,RebarAvailable,arreglo_t1,...
                            arreglo_t2,pu_steel_footings,wac);
                        
        bestCost_elem=bestCost_elem+bestCost_axis;
        dispositionRebarFinal1=dispositionRebar;
        arreglo_t1B=arreglo_t1;
        arreglo_t2B=arreglo_t2;
    else
        acmin=pmin*be*d;
        areaBar=ac+acmin;
        AcBar=[AcBar,areaBar];
        
        list_ef_footings(2)=bestEf;    % section Le
        list_mr_footings(2)=bestMr;
        
        barDispositionFootings=[];

        arrangement_bar_footings=[];

        nbars_footings(1,2*eje-1)=length(arreglo_t1);
        nbars_footings(1,2*eje)=length(arreglo_t2);
        bestCost_axis = EvaluateCostISRFoot(le,rec,ac,acmin,...
            pu_steel_footings,wac);
                    
        bestCost_elem=bestCost_elem+bestCost_axis;
    end
    
end
if PlotRebarDesign==1
    ReinforcedSectionsFooting(hmodif,be,le,dispositionRebarFinal1,...
    dispositionRebarFinal2,arreglo_t1B,arreglo_t2B,arreglo_t1L,arreglo_t2L);
end