function [totalWeightStruc,wsteelColsTotal,pacColsElem,Mp,dimensions,...
    unitWeightElem,wsteelConcBeamsElem,wsteelConcColsElem,...
    wsteelConcFootingsElem,hefootings,dimFoot,totalCostStruc,inertiaElem,...
    wsteelStructure]=DesignRCPlaneFrameBCI(puBeams,puCols,lenElem,fcElem,...
    inertiaElem,qadm,FSfootings,nodesSupportColumns,puSteelFootings,...
    dimensions,fcbeams,fccols,fcfootings,areaElem,colsSymAsymISR,...
    rebarAvailable,conditionCracking,duct,elem_cols,elem_beams,recCols,...
    load_conditions_beams,load_conditions_columns,reactions,shearbeam,...
    coordBaseCols,coordEndBeams,coordBaseFooting,directionData)
%% Documentation:
%------------------------------------------------------------------------
% Syntax:
% [totalWeightStruc,wsteelColsTotal,pacColsElem,Mp,dimensions,...
%  unitWeightElem,wsteelConcBeamsElem,wsteelConcColsElem,...
%  wsteelConcFootingsElem,hefootings,dimFoot,totalCostStruc,inertiaElem,...
%  wsteelStructure]=DesignRCPlaneFrameBCI(puBeams,puCols,lenElem,fcElem,...
%  inertiaElem,qadm,FSfootings,nodesSupportColumns,puSteelFootings,...
%  dimensions,fcbeams,fccols,fcfootings,areaElem,colsSymAsymISR,...
%  rebarAvailable,conditionCracking,duct,elem_cols,elem_beams,recCols,...
%  load_conditions_beams,load_conditions_columns,reactions,shearbeam,...
%  coordBaseCols,coordEndBeams,coordBaseFooting,directionData)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To optimally design the reinforcement of the elements of a 
% reinforced concrete plane frame composed of rectangular beams, 
% rectangular columns and rectangular isolated footings.
% 
% OUTPUT: totalWeightStruc:             is the total weight of the structural
%                                       frame, considering the reinforcing 
%                                       steel and the concrete volume
%
%         wsteelColsTotal:              is the sum of reinforcing steel 
%                                       weight of each of the columns 
%                                       composing the structural frame
%
%         pacColsElem:                  is a vector containing the percentage
%                                       steel area of each column. 
%                                       Size = [ncolumns,1]
%
%         Mp:                           are the resisting bending moments 
%                                       at the ends of each element composing
%                                       the structural frame (Plastic Moments)
%
%         dimensions:                   are the new modified cross-section
%                                       dimensions of the elements
%
%         unitWeightElem:               Is the array containing the unit-
%                                       self-weight of each structural
%                                       element considering both the
%                                       concrete and steel reinforcement
%
%         wsteelConcBeamsElem:          is the vector containing the total
%                                       weight of each of the beam elements,
%                                       considering both the concrete volume
%                                       and the steel reinforcement
%
%         wsteelConcColsElem:           is the vector containing the total
%                                       weight of each of the column elements,
%                                       considering both the concrete volume
%                                       and steel reinforcement
%
%         wsteelConcFootingsElem:       is the vector containing the total
%                                       weight of each of the footing 
%                                       elements, considering both the
%                                       concrete volume and steel 
%                                       reinforcement
%
%         hefootings:                   is the vector containing the final
%                                       designed height dimensions of the
%                                       isolated footings
%
%         dimFoot:                      is the vector containing the 
%                                       transversal cross-section dimensions
%                                       of all footings. Size=[nfootings,2] 
%                                       in format [Be,Le]
%
%         totalCostStruc:               Is the total construction cost of
%                                       the structural frame considering the
%                                       reinforcing steel design and
%                                       concrete volumes
%
%         inertiaElem:                  are the final cross-section inertia
%                                       momentums of each structural element
%                                       after applying a cracked or 
%                                       non-cracked cross-section mechanism
%
%         wsteelStructure:              is the total weight of steel
%                                       reinforcement of the structural
%                                       frame
%
% INPUT:  puBeams:                      is the unit construction cost of
%                                       steel reinforcement assembly for a
%                                       beam element, as an average of all
%                                       assembly performances for each type
%                                       of rebar commercially available
%                                       (considering that as many as 6
%                                       different types of rebar may be
%                                       placed in a beam element)
%
%         puCols:                       are the unit construction cost data
%                                       for the reinforcing bar assembly.
%                                       For symmetric reinforcement format
%                                       is default: 
%    _________________________________________________________________
%    pu_col=[PU{#4}, PU{#5}, PU{#6}, PU{#8}, PU{#9}, PU{#10}, PU{#12}]
%    -----------------------------------------------------------------
%                                       and for asymmetric reinforcement as
%                                       an array of size: [2,7] by default
%                                       for asymmetric reinforcement, for
%                                       which the first row corresponds to 
%                                       unit-cost values for each rebar 
%                                       type in case only one type of rebar
%                                       results as an optimal design, and
%                                       the second row consisting of only
%                                       one unit cost in case more than one
%                                       different type of rebar results as 
%                                       an optimal design (assuming an 
%                                       average unit-cost of all types of 
%                                       rebars)
%
%         lenElem:                      is an array containing the length
%                                       of each structural element composing
%                                       a structural frame
%
%         fcElem:                       is the vector containing the f'c 
%                                       used for each of the elements 
%                                       composing the structural frame
%
%         inertiaElem:                  is the vector containing the cross-
%                                       section inertia momentum of each of
%                                       the elements composing the 
%                                       structural frame
%
%         qadm:                         is the max admissible bearing load 
%                                       of the soil supporting the footings
%                                       of the structural frame. Is used for
%                                       the design of the footings
%
%         FSfootings:                   is the Safety Factor used for the 
%                                       design of the footings
%
%         nodesSupportColumns:          is an array containing the nodes in
%                                       which each of the base columns are 
%                                       supported in contact with their
%                                       respective isolated footing. 
%                                       Size = [2,nfootings] in format:
%
%                                        [node;
%                                       columnID]
%
%         puSteelFootings:              is the unit construction cost of
%                                       steel reinforcement assembly for an 
%                                       isolated footing element, as an 
%                                       average of all assembly performances
%                                       for each type of rebar commercially
%                                       available (considering that as many
%                                       as 4 different types of rebar may be
%                                       placed in an isolated footing element)
%
%         dimensions:                   is the array containing the cross-
%                                       section dimensions of all elements
%                                       composing the structural frame. 
%                                       Size = [nbars,2] in format [be,he]
%
%         fcbeams,fccols,fcfootings:    is the f'c used for beam, columns 
%                                       and isolated footings elements,
%                                       respectively 
%
%         areaElem:                     is the vector containing the cross-
%                                       section area of all of the structural
%                                       elements composing the structural
%                                       frame
%
%         ForcesDOFseismic:             is the list of DOF on which the 
%                                       lateral equivalent base shear forces
%                                       are acting
%
%         colsSymAsymISR:               is the parameter that indicates what
%                                       type of reinforcement design is to
%                                       be carried out (mainly for columns,
%                                       although it also defines the type 
%                                       of design for the other element 
%                                       types). Options are ''ISR'',
%                                       ''Symmetric'',''Asymmetric''
%
%         conditionCracking:            is the parameter that indicates what
%                                       type of cracking mechanism is 
%                                       considered for the computation of 
%                                       the inertia momentums of each of the 
%                                       elements' cross-section. Options are:
%                                       ''Cracked'' and ''Non-cracked''
%
%         elem_cols,elem_beams:         are the vector containing the 
%                                       elements' IDs or numbers 
%                                       corresponding to the type of 
%                                       structural element. Only for columns 
%                                       and beams, respectively
%
%         recxyCols:                    is a vector containing the concrete
%                                       cover for both axis directions of a
%                                       rectangular column cross-section as
%                                       coverX,coverY
%
%         directionData:                is the disc route direction to save 
%                                       the design data (if required). 
%                                       Two options are available:
%
%                                       - an empty vector [] if it is not
%                                         required to export the data
%                                       - the route direction as a string
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
%% ----------------------------------------------------------------------
%                  OPTIMAL DESIGN OF REBAR ON ELEMENTS
% -----------------------------------------------------------------------

%% DESIGN OF BEAMS 
nbeams=length(elem_beams);
                
h_rec_beams=4; % concrete cover for beams
b_rec_beams=4; % (cm)                 
fy=4200; % yield stress of steel reinforcement
wac=7.8e-3; % unit volume weight of the reinforcement (Kg/cm3)
wco=2.4e-3; % unit volume weight of the concrete (Kg/cm3)
graphConvergenceBeamPlot=0;
rebarBeamDesignPlots=0;

sectionRestrictions=0;
disposition_rebar_beams3sec=[];
typeRebarsBeams=[];
for i=1:nbeams
    nelem=elem_beams(i);
    b=dimensions(nelem,1); %cm
    h=dimensions(nelem,2); %cm
    span=lenElem(nelem); %cm
    fc=fcElem(nelem);

    [sepbarsRestric,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
    barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
    barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
    barArrangementIzqComp,minAreabar3sec,Efelem3secISR,bestCostVar,efrebar3sec,...
    minAreaRebarAverage,Mr3section]=beamsISR(puBeams,span,wac,b,h,h_rec_beams,...
    rebarAvailable,fc,fy,load_conditions_beams(i,:),colsSymAsymISR,duct,b_rec_beams,...
    rebarBeamDesignPlots,graphConvergenceBeamPlot);

    if h/b>3 || span/h<=4 || sum(sepbarsRestric)~=0
        sectionRestrictions=1; % to know if the geometry 
                               % restrictions comply 
    end

    % To collect and save design results:
    MrbeamsCollection(i,:)=Mr3section;
    inertiaModifbeams(i,:)=inertia_modif;
                             
    dim_beams_collection(i,:)=[b h span h_rec_beams];
             
    disposition_rebar_beams3sec=[disposition_rebar_beams3sec;
                                dispositionBar_Izq;
                                 dispositionBar_Center;
                                 dispositionBar_Der];
                         
    nbarbeamsCollection(i,:)=[length(barArrangementIzqTens)...
                             length(barArrangementIzqComp)...
                             length(barArrangementCentralTens)...
                             length(barArrangementCentralComp)...
                             length(barArrangementDerTens)...
                             length(barArrangementDerComp)];  
    
    typeRebarsBeams=[typeRebarsBeams;
                    barArrangementIzqTens;
                    barArrangementIzqComp;
                    barArrangementCentralTens;
                    barArrangementCentralComp;
                    barArrangementDerTens;
                    barArrangementDerComp];
                
    % Update dimensions of beams in the global model (in case they were
    % modified during the design process):
    dimensions(nelem,2)=h; %cm
    %% SHEAR DESIGN
    rho=(sum(minAreabar3sec)/length(minAreabar3sec))/(b*h);
    
    [s1,s2,s3,d1,d3]=shearDesignBeams(span,b,h,h_rec_beams,rho,fc,fy,...
                     shearbeam(i,:));

    shear_beam_design_collec(i,:)=[s1,s2,s3,d1,d3];  
    
    % Additional relevant design results:
    if colsSymAsymISR~="ISR"
        efrebarbeams(i,:)=efrebar3sec;
        efbeamsAverage(i,1)=sum(efrebarbeams(i,:))/3;
    else
        efrebarbeams(i,:)=Efelem3secISR;
        efbeamsAverage(i,1)=sum(efrebarbeams(i,:))/3;
    end
    costSteelbeams(i,1)=bestCostVar;
    areaSteelbeams3sec(i,:)=minAreabar3sec;
    AreaRebarAverageBeams(i,1)=minAreaRebarAverage;

    % To compute the resistance at each end of the beam elements
    Mp(elem_beams(i),1)=MrbeamsCollection(i,1);
    Mp(elem_beams(i),2)=MrbeamsCollection(i,3);

    % To modify the initial given "inertia" vector with the new modified
    % ones 
    inertiaElem(elem_beams(i))=sum(inertiaModifbeams(i,:))/3;
end

% Exporting results - if required
if isempty(directionData)==0
    ExportResultsBeam(directionData,dim_beams_collection,coordEndBeams,...
    disposition_rebar_beams3sec,nbarbeamsCollection,...
    typeRebarsBeams,colsSymAsymISR,shear_beam_design_collec);

end
%% DESIGN OF COLUMNS
ncols=length(elem_cols);

optimaColConvPlot=0;
plotISRColResults=0;
plotRebarDesign=0;


nlf=length(load_conditions_columns(:,1))/ncols; % This variable is to
                                                % to know how many loads
                                                % take place for each 
                                                % column.
    
bestdisposicionRebarCollec=[];
typeRebarsColsCollection=[];
for i=1:ncols
    nelem=elem_cols(i);
    b=dimensions(nelem,1); %cm
    h=dimensions(nelem,2); %cm
    
    height=lenElem(nelem); %cm
    fc=fcElem(nelem); %kg/cm2
    
    loadCol=load_conditions_columns(nlf*i-(nlf-1):nlf*i,:);
    
    [Inertia_xy_modif,b,h,bestArrangement,bestdisposicionRebar,costCol,...
    AcSecCols,EfsecCol,Mr_col]=isrColumnsSymAsym(puCols,height,b,h,recCols,...
    fy,fc,loadCol,colsSymAsymISR,conditionCracking,duct,wac,rebarAvailable,...
    optimaColConvPlot,plotISRColResults,plotRebarDesign);
    

    if b/h<0.4 || height/b>15
        sectionRestrictions=1;
    end
    MrColumnsCollection(i,:)=Mr_col;
    inertia_x_modif_cols(i,1)=Inertia_xy_modif(1);
    
    dimensions(nelem,2)=h; % in case hte height dimension was modified 
    
    % -------------------------------------------------------------------
    % Collecting design results of all columns 
    % -------------------------------------------------------------------
    dimColumnsCollection(i,:)=[b h height recCols(1)];
                                    
    nbarColumnsCollection(i,1)=length(bestArrangement);
    
    typeRebarsColsCollection=[typeRebarsColsCollection;
                            bestArrangement];
                
    bestdisposicionRebarCollec=[bestdisposicionRebarCollec;
                                bestdisposicionRebar];
    
    % -------------------------------------------------------------------
    
    steelAreacols(i,1)=AcSecCols;
    cost_steel_cols(i,1)=costCol;
    ef_cols(i,1)=EfsecCol;

    inertiaElem(elem_cols(i))=inertia_x_modif_cols(i,1);
    Mp(elem_cols(i),1)=MrColumnsCollection(i,1);
    Mp(elem_cols(i),2)=MrColumnsCollection(i,1);
end

% Exporting results - if required
if isempty(directionData)==0
    ExportResultsColumn(directionData,dimColumnsCollection,...
    bestdisposicionRebarCollec,nbarColumnsCollection,...
    typeRebarsColsCollection,colsSymAsymISR,coordBaseCols);

end
%% ----------------------------------------------------------------------
% DESIGN OF FOOTINGS 
% -----------------------------------------------------------------------

qu=qadm/FSfootings;
recFoot=5;

nfootings=length(nodesSupportColumns(1,:));
hefoot=30; % 30 cm for all footings, initially
typeFoot=1; % isolated standard footings (with column at the center)
emin=2; %min load eccentricity in the out-of-plane direction (cm)

optimConvFootings=0;
PlotRebarFootingDesign=0;

collectionArrangementRebarFootings=[];
collectionDispositionRebarFootings=[];
for j=1:nfootings
    napoyo=nodesSupportColumns(1,j);
    
    % Extract design loads from the structural analysis data - reactions
    pu=-reactions(3*nodesSupportColumns(1,j)-1);
    mux=-reactions(3*nodesSupportColumns(1,j));
    
    col_apoyo=nodesSupportColumns(2,j);
    dimCol=dimensions(col_apoyo,:);
    
    % Dtermine initial transversal dimensions of the footing
    [be,le,contact_pressures(j)]=designDimFootings(pu,qu,dimCol,hefoot,...
        recFoot,typeFoot);
    
    dimFoot(j,2)=le;
    dimFoot(j,1)=be;
    
    load_conditions_footings(j,1)=napoyo;
    load_conditions_footings(j,2)=pu;
    load_conditions_footings(j,3)=mux;
    load_conditions_footings(j,4)=pu*emin;
    
    % Optimal design of reinforcement in the Footing
    [hmodif,mu_axis,bestDispositionFootings,arrangement_bar_footings,...
    nbars_footings,AcBar,bestCost_elem,Ef_elemBar,list_mr_footings]=...
    isrFootings(puSteelFootings,hefoot,be,le,recFoot,fc,fy,...
    load_conditions_footings(j,:),dimCol,rebarAvailable,colsSymAsymISR,...
    duct,optimConvFootings,PlotRebarFootingDesign,typeFoot,10,wac);

    hefootings(j,1)=hmodif;
    
    % Collecting desing results for exportation (if required)
    collectionArrangementRebarFootings=[collectionArrangementRebarFootings;
                                        arrangement_bar_footings];
    
    collectionDispositionRebarFootings=[collectionDispositionRebarFootings;
                                        bestDispositionFootings];
                                    
              
    
    nbarsFootingsCollection(j,:)=nbars_footings;
    
    dimensionFootingCollection(j,:)=[be le hmodif recFoot];
    
    areaElemFootingsRebar(j,:)=AcBar;
    efElemFootingsRebar(j,:)=Ef_elemBar;
    costElemFootings(j,1)=bestCost_elem;
end

% Exporting design results (if required)
if isempty(directionData)==0
    ExportResultsIsolFootings(directionData,collectionDispositionRebarFootings,...
    dimensionFootingCollection,nbarsFootingsCollection,collectionArrangementRebarFootings,...
    coordBaseFooting,colsSymAsymISR)
end

%% WEIGHT OF STRUCTURE 
for i=1:nbeams
    areaSteelbeams(i,1)=sum(areaSteelbeams3sec(i,:));
end

[wsteelColsTotal,pacColsElem,wsteelConcBeamsElem,wsteelConcColsElem,...
wsteelConcFootingsElem,vol_beams,vol_cols,vol_footings,wsteelStructure,...
wconc_structure,wbeams,wcols,wfootings,totalWeightStruc]=WeightStruc...
(elem_cols,elem_beams,lenElem,areaElem,areaSteelbeams,areaElemFootingsRebar,...
hefootings,nbeams,ncols,steelAreacols,nfootings,dimFoot,wac,wco);

% The unit weight of each element is updated so that the contribution of
% the reinforcing steel is considered:

for i=1:nbeams
    unitWeightElem(elem_beams(i),2)=wsteelConcBeamsElem(i)/...
                        (areaElem(elem_beams(i))*lenElem(elem_beams(i)));
    unitWeightElem(elem_beams(i),1)=elem_beams(i);
end

for i=1:ncols
    unitWeightElem(elem_cols(i),2)=wsteelConcColsElem(i)/...
                            (areaElem(elem_cols(i))*lenElem(elem_cols(i)));
    unitWeightElem(elem_cols(i),1)=elem_cols(i);

end

%% COST OF STRUCTURE
             %f'c % Bombed %Direct Shot (per unit volume) $/cm3
cost_concrete=[100 2281.22e-6 2266.98e-6;
                150 2401.22e-6 2390.88e-6;
                200 2532.14e-6 2525.28e-6;
                250 2777.96e-6 2845.00e-6;
                300 2939.12e-6 3010.90e-6;
                350 3111.50e-6 3188.50e-6;
                400 3298.16e-6 3380.50e-6;
                450 3499.10e-6 3587.35e-6];
            
totalCostStruc=CostStruc(costSteelbeams,cost_steel_cols,costElemFootings,...
fcbeams,fccols,fcfootings,vol_beams,vol_cols,vol_footings,cost_concrete);

