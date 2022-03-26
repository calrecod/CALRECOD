function [totalWeightStruc,wsteelColsTotal,pacColsElem,sectionRestrictions,...
    Mp,dimensions,displacements,unitWeightElem,wsteelConcBeamsElem,...
    wsteelConcColsElem,wsteelConcFootingsElem,hefootings,dimFoot,...
    totalCostStruc,inertiaElem,wsteelStructure]=DesignRCPlaneFrameBCIF...
    (nnodes,bc,ni,nf,Edof,coordxy,puBeams,type_elem,puCols,nbars,np,...
    lenElem,coordBaseCols,fcElem,inertiaElem,qadm,FSfootings,nodesSupportColumns,...
    puSteelFootings,dimensions,Eelem,fcbeams,fccols,fcfootings,fglobal,...
    qbary,areaElem,ForcesDOFseismic,floorElem,colsSymAsymISR,...
    conditionCracking,duct,elem_cols,elem_beams,recxyCols,directionData,...
    plotAnalysisResults)

%------------------------------------------------------------------------
% Syntax:
% [totalWeightStruc,wsteelColsTotal,pacColsElem,sectionRestrictions,...
%   Mp,dimensions,displacementsRightLeft,unitWeightElem,wsteelConcBeamsElem,...
%   wsteelConcColsElem,wsteelConcFootingsElem,hefootings,dimFoot,...
%   totalCostStruc,inertiaElem,wsteelStructure]=DesignRCPlaneFrameBCIF...
%   (nnodes,bc,ni,nf,Edof,coordxy,puBeams,type_elem,puCols,nbars,np,...
%   lenElem,coordBaseCols,fcElem,inertiaElem,qadm,FSfootings,nodesSupportColumns,...
%   puSteelFootings,dimensions,Eelem,fcbeams,fccols,fcfootings,fglobal,...
%   qbary,areaElem,ForcesDOFseismic,floorElem,colsSymAsymISR,...
%   conditionCracking,duct,elem_cols,elem_beams,recxyCols,directionData,...
%   plotAnalysisResults)
%
%------------------------------------------------------------------------
% PURPOSE: To design the reinforcement of the elements of a plane frame
% concrete structure composed of rectangular beams, rectangular columns 
% and rectangular isolated footings.
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
%         sectionRestrictions:          is a parameter that indicates weather
%                                       or not the dimension constraints of
%                                       both beam and column elements are
%                                       being complied. If the parameter's
%                                       value is equal to 1, then the 
%                                       dimension restrictions are not being
%                                       complied, if its value is equal to 0
%                                       all dimension restrictions are being
%                                       complied
%
%         Mp:                           are the resisting bending moments 
%                                       at the ends of each element composing
%                                       the structural frame (Plastic Moments)
%
%         dimensions:                   are the new modified cross-section
%                                       dimensions of the elements
%
%         displacementsRightLeft:       is an array containing the node
%                                       displacements of the structural
%                                       frame extracted from both structural
%                                       static linear analysis with lateral
%                                       seismic loads to the right and to
%                                       the left. Size = [NDOF,2] in format: 
%                                       [diplace{right},diplace{left}]
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
% INPUT:  nnodes:                       is the number of nodes of the 
%                                       structural frame
%
%         bc:                           is the boundary condition vector of 
%                                       the restricted DOF with their pre-
%                                       established displacement condition. 
%                                       Size = [DOF,2] in format: 
%                                       [nDOF,displacement]
%
%         ni,nf:                        are the list of initial and final 
%                                       nodes ID's, respectively, for each
%                                       of the elements composing a 
%                                       structural frame in the order pre-
%                                       established
%
%         Edof:                         is the topology array.
%                                       Size=[nbars,7] in format:  
%               ____________________________________________
%               [no{bar},DOF{initial-node},DOF{final-node}]
%               --------------------------------------------
%
%         coordxy:                      is an array containing the node
%                                       coordinates. Size = [n{nodes},2] in
%                                       format: [x,y]
%
%         puBeams:                      is the unit construction cost of
%                                       steel reinforcement assembly for a
%                                       beam element, as an average of all
%                                       assembly performances for each type
%                                       of rebar commercially available
%                                       (considering that as many as 6
%                                       different types of rebar may be
%                                       placed in a beam element)
%
%         type_elem:                    is an array containing the element
%                                       ID labels of each of the elements
%                                       composing a structural frame: either
%                                       ''Beam'' or ''Col''. 
%                                       Size = [nbars,2] in format: 
%                                       [nbar,''label'']
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
%         nbars:                        is the number of bars or elements
%                                       composing a structural frame
%
%         np:                           are the number of points of analysis
%                                       for the computation of the mechanic
%                                       elements distribution of an element
%
%         lenElem:                      is an array containing the length
%                                       of each structural element composing
%                                       a structural frame
%
%         coordBaseCols:                are the coordinates of each of the
%                                       columns' base point. 
%                                       Size = [ncolumns,2] in format: 
%                                       [x,y]
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
%         Eelem:                        is the yield stress of each of the
%                                       structural elements' material (is a
%                                       function of the $f'c used for each
%                                       element)
%
%         fcbeams,fccols,fcfootings:    is the f'c used for beam, columns 
%                                       and isolated footings elements,
%                                       respectively 
%
%         fglobal:                      is the vector containing the 
%                                       punctual node forces applied to the 
%                                       structural frame. Size = [NDOF,1]
%
%         qbary:                        is the array containing the value
%                                       of the distributed load acting
%                                       downwards over each of the elements 
%                                       (applies specially for beams)
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
%         floorElem:                    is an array containing the elements'
%                                       IDs that correspond to each of the 
%                                       floors of the structural element. 
%                                       Size = [nfloors,nmaxElements+1]. 
%                                       It is used to compute the weight, 
%                                       mass and stiffness of each floor. 
%                                       Format: 
%                   ______________________________
%                   [nfloor,elem_{1},...,elem_{n}]
%                   ------------------------------
%                                       In case one floor has less number 
%                                       of elements than the others, the 
%                                       array should be completed by zeros
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
%                                       the design data
%
%         plotAnalysisResults:          is the parameters that indicates if
%                                       the linear static structural analysis
%                                       results are required or not. Options 
%                                       are: (1) they are required, (2) they
%                                       are not required
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

Mp=zeros(nbars,2);

nbeams=length(elem_beams);
ncols=length(elem_cols);

% TWO ANALYSIS AS EXECUTED IN CASE LATERAL EQUIVALENT SEISMIC LOADS
% ARE APPLIED TO THE FRAME_____________________________________________
if isempty(ForcesDOFseismic)==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% analysis taking all loads in consideration (F-der)............
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [displacements_right,r_right,Ex_right,Ey_right,es_bars_normal_right,...
    es_bars_shear_right,es_bars_moment_right]=PlaneFrameStaticLinearAnalysis(nnodes,...
    nbars,Eelem,areaElem,inertiaElem,bc,fglobal,ni,nf,qbary,Edof,np,coordxy,plotAnalysisResults);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% analysis taking all loads in consideration (F-izq)................
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fglobal_izq=fglobal;
    fglobal_izq(ForcesDOFseismic)=-fglobal(ForcesDOFseismic);


    [displacements_left,r_left,Ex_left,Ey_left,es_bars_normal_left,...
    es_bars_shear_left,es_bars_moment_left]=PlaneFrameStaticLinearAnalysis(nnodes,...
    nbars,Eelem,areaElem,inertiaElem,bc,fglobal_izq,ni,nf,qbary,Edof,np,coordxy,plotAnalysisResults);

    displacements=[displacements_right,displacements_left];
else
    [displacements_right,r_right,Ex_right,Ey_right,es_bars_normal_right,...
    es_bars_shear_right,es_bars_moment_right]=PlaneFrameStaticLinearAnalysis(nnodes,...
    nbars,Eelem,areaElem,inertiaElem,bc,fglobal,ni,nf,qbary,Edof,np,coordxy,plotAnalysisResults);

    displacements=[displacements_right];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% DESIGN OF ELEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mechanic elements  __________________________________________________
es_bars_normal_elems=es_bars_normal_right;
es_bars_shear_elems=es_bars_shear_right;
es_bars_moment_elems=es_bars_moment_right;

if isempty(ForcesDOFseismic)==0
    % To determine the max mechanic element values on each element for each
    % analysis point (compare absoulute values) ___________________________
    for j=1:nbars
        for i=1:np
            if abs(es_bars_normal_elems(i,j))<abs(es_bars_normal_left(i,j))
                es_bars_normal_elems(i,j)=es_bars_normal_left(i,j);

            end
            if abs(es_bars_shear_elems(i,j))<abs(es_bars_shear_left(i,j))
                es_bars_shear_elems(i,j)=es_bars_shear_left(i,j);

            end
            if abs(es_bars_moment_elems(i,j))<abs(es_bars_moment_left(i,j))
                es_bars_moment_elems(i,j)=es_bars_moment_left(i,j);

            end
        end
    end
end
unitWeightElem=zeros(nbars,2);

% DESIGN OF BEAMS _____________________________________________________
% ______________________________________________________________________

es_bars_normal_elems=es_bars_normal_elems*0.001; % ton
es_bars_moment_elems=es_bars_moment_elems*0.00001; %ton-m
     
% To determine the load conditions at each of the three critical design
% section of the beam __________________________________________________

load_conditions_beams=zeros(nbeams,4);
for i=1:nbeams
    nelem=elem_beams(i);
    load_conditions_beams(i,1)=es_bars_normal_elems(1,nelem);
    load_conditions_beams(i,2)=es_bars_moment_elems(1,nelem);
    load_conditions_beams(i,3)=max(abs(es_bars_moment_elems(3:5-1,nelem)));
    load_conditions_beams(i,4)=es_bars_moment_elems(7,nelem);
end

inertiaModifbeams=zeros(nbeams,3);
efrebarbeams=zeros(nbeams,3);
costSteelbeams=zeros(nbeams,1);
areaSteelbeams3sec=zeros(nbeams,6);
efbeamsAverage=zeros(nbeams,1);
AreaRebarAverageBeams=zeros(nbeams,1);

collectionDispositionRebarBeams3sec=[];
collectionArrangebarbeams=[];
dim_beams_collection=zeros(nbeams,9);
nbarbeamsCollection=zeros(nbeams,6);
                
h_rec_sections=[5 5 3 3 5 5]; % [rec_left_up, rec_left_low, rec_mid_up,
                              % rec_mid_low, rec_right_up, rec_right_low]
b_rec=4; % (cm)                 
fy=4200;

% To export design results ___________________________________________
nombre_archivo='nrebar_beams_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_beams_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_beams.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_beams3sec.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');
%______________________________________________________________________

graphConvergenceBeamPlot=0;
rebarBeamDesignPlots=0;

MrbeamsCollection=zeros(nbeams,3);
sectionRestrictions=0;
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
    minAreaRebarAverage,Mr3section]=beamsISR(puBeams,span,b,h,h_rec_sections,...
    fc,fy,load_conditions_beams,colsSymAsymISR,duct,b_rec,rebarBeamDesignPlots,...
    graphConvergenceBeamPlot);

    if h/b>3 || span/h<=4 || sum(sepbarsRestric)~=0
        sectionRestrictions=1; % to know if the restrictions comply 
    end

    % To collect and save design results _______________________________
    MrbeamsCollection(i,:)=Mr3section;
    inertiaModifbeams(i,:)=inertia_modif;
    disposition_rebar_beams3sec=[dispositionBar_Izq;
                                 dispositionBar_Center;
                                 dispositionBar_Der];
                             
    collectionDispositionRebarBeams3sec=[collectionDispositionRebarBeams3sec;
                                            disposition_rebar_beams3sec];
    
    dim_beams_collection(i,:)=[b h span h_rec_sections];
    
    arrangemetbarbeams=[barArrangementIzqTens;
                        barArrangementIzqComp;
                        barArrangementCentralTens;
                        barArrangementCentralComp;
                        barArrangementDerTens;
                        barArrangementDerComp];
                    
    collectionArrangebarbeams=[collectionArrangebarbeams;
                                    arrangemetbarbeams];
                                
    nbarbeamsCollection(i,:)=[length(barArrangementIzqTens)...
        length(barArrangementIzqComp) length(barArrangementCentralTens)...
        length(barArrangementCentralComp) length(barArrangementDerTens)...
        length(barArrangementDerComp)];  
    
    if colsSymAsymISR~="ISR"
        efrebarbeams(i,:)=efrebar3sec;
        efbeamsAverage(i)=sum(efrebarbeams(i,:))/3;
    else
        efrebarbeams(i,:)=Efelem3secISR;
        efbeamsAverage(i)=sum(efrebarbeams(i,:))/3;
    end
    costSteelbeams(i)=bestCostVar;
    areaSteelbeams3sec(i,:)=minAreabar3sec;
    AreaRebarAverageBeams(i)=minAreaRebarAverage;

    % To compute the resistance at each end of the beam elements _______
    Mp(elem_beams(i),1)=MrbeamsCollection(i,1);
    Mp(elem_beams(i),2)=MrbeamsCollection(i,3);

    % To modify the initial given "inertia" vector with the new modified
    % ones _____________________________________________________________
    inertiaElem(elem_beams(i))=sum(inertiaModifbeams(i,:))/3;
    
    % To write .txt files for the exportation of results _____________
    if colsSymAsymISR~="ISR"
        fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',...
                        dim_beams_collection(i,:));
    
    
        for j=1:length(disposition_rebar_beams3sec(:,1))
            fprintf(fileid_04,'%.2f %.2f\n',disposition_rebar_beams3sec(j,:));
        end
        
        fprintf(fileid_01,'%d %d %d %d %d %d\n',nbarbeamsCollection(i,:));
                  
        fprintf(fileid_03,'%d\n',arrangemetbarbeams(:,:));
    else
        fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',...
                        dim_beams_collection(i,:));
        
        fprintf(fileid_04,'%.2f %.2f\n',[0,0]);
        
            
        fprintf(fileid_01,'%d %d %d %d %d %d\n',[0,0,0,0,0,0]);

        fprintf(fileid_03,'%d\n',[0]);
    end
end
fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

% DESIGN OF COLUMNS ____________________________________________________
% ______________________________________________________________________

emin=0.02; %minimum excentricity (m)
load_conditions_columns=zeros(ncols,4);
for i=1:ncols
    nelem=elem_cols(i);
    axial=es_bars_normal_elems(1,nelem);

    Mux=max(es_bars_moment_elems(:,elem_cols(i)))*0.00001;
    Mminy=abs(axial*emin);
    load_conditions_columns(i,2)=axial;
        
    load_conditions_columns(i,3)=Mux; % Mux (in-plane)
    load_conditions_columns(i,4)=Mminy; % Muy (out-of-plane)
end

inertia_x_modif_cols=zeros(ncols,1);

% To export design results __________________________________________
nombre_archivo='nrebar_cols_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_cols_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_cols.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_cols.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');
% ____________________________________________________________________
optimaColConvPlot=0;
plotISRColResults=0;
plotRebarDesign=0;

steelAreacols=zeros(ncols,1);
ef_cols=zeros(ncols,1);
cost_steel_cols=zeros(ncols,1);

collectionArragnementCols=[];
collectionDispositionRebarCols=[];
nbarColumnsCollection=zeros(ncols,1);
dimColumnsCollection=zeros(ncols,5);
MrColumnsCollection=zeros(ncols,2);
for i=1:ncols
    nelem=elem_cols(i);
    b=dimensions(nelem,1); %cm
    h=dimensions(nelem,2); %cmteclear
    rec=recxyCols;
    
    height=lenElem(nelem); %cm
    fc=fcElem(nelem); %kg/cm2
    
    [Inertia_xy_modif,b,h,bestArrangement,bestdisposicionRebar,costCol,...
    AcSecCols,EfsecCol,Mr_col]=isrColumnsSymAsym(puCols,height,b,h,rec,...
    fy,fc,load_conditions_columns,colsSymAsymISR,conditionCracking,duct,...
    optimaColConvPlot,plotISRColResults,plotRebarDesign);
    
    if b/h<0.4 || height/b>15
        sectionRestrictions=1;
    end
    
    MrColumnsCollection(i,:)=Mr_col;
    inertia_x_modif_cols(i)=Inertia_xy_modif(1);
    
    dimensions(nelem,2)=h; %cmteclear
    
    dimColumnsCollection(i,:)=[b h height rec];
    
    collectionArragnementCols=[collectionArragnementCols;
                                bestArrangement];
                                    
    
    collectionDispositionRebarCols=[collectionDispositionRebarCols;
                                         bestdisposicionRebar];
    
    
    nbarColumnsCollection(i)=length(bestArrangement);
    
    steelAreacols(i)=AcSecCols;
    cost_steel_cols(i)=costCol;
    ef_cols(i)=EfsecCol;

    inertiaElem(elem_cols(i))=inertia_x_modif_cols(i);
    Mp(elem_cols(i),1)=MrColumnsCollection(i,1);
    Mp(elem_cols(i),2)=MrColumnsCollection(i,1);
    
    % Writing the .txt file for the exportation of results______________
    if colsSymAsymISR~="ISR"
        fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f\n',dimColumnsCollection(i,:));
                               
        fprintf(fileid_03,'%d\n',bestArrangement);
    
        for j=1:length(bestdisposicionRebar(:,1))
            fprintf(fileid_04,'%.2f %.2f\n',bestdisposicionRebar(j,:));
        end
    
        fprintf(fileid_01,'%d\n',nbarColumnsCollection(i));
    else
        fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f\n',dimColumnsCollection(i,:));
                               
        fprintf(fileid_03,'%d\n',bestArrangement);
    
    end
    % _________________________________________________________________
end
fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

% To adjust cross-section dimensions of elements _______________________

% Lower columns are to have smaller o equal cross-section dimensions than
% the upper ones .......................................................

nfloors=length(floorElem(:,1));
for i=1:2*nfloors
    [dimensions]=AdjustDimElemFrames(nbars,type_elem,nnodes,dimensions,...
        coordBaseCols,ni,nf);
end 

%_______________________________________________________________________

% DESIGN OF FOOTINGS ___________________________________________________
% ______________________________________________________________________
if isempty(ForcesDOFseismic)==0
    for j=1:length(r_right)
        if abs(r_right(j))<abs(r_left(j))
            r_right(j)=r_left(j);
        end

    end
end
reactions=r_right;

qu=qadm*FSfootings;
rec=5;

nzapatas=length(nodesSupportColumns(1,:));
hefootings=zeros(nzapatas,1)+30; % 30 cm for all footings initially

emin=0.02; %m
load_conditions_footings=zeros(nzapatas,4);
contact_pressures=zeros(nzapatas,1);
for j=1:nzapatas
    napoyo=nodesSupportColumns(1,j);
    
    pu=reactions(3*nodesSupportColumns(1,j)-1);
    muy=abs(reactions(3*nodesSupportColumns(1,j)));
    
    col_apoyo=nodesSupportColumns(2,j);
    dimCol=dimensions(col_apoyo,:);
    
    lebe=sqrt(pu/qu);
    
    hfooting=hefootings(j);
    be=lebe;
    le=lebe;
    
    if (be-dimCol(1,2))*0.5<((hfooting-rec)*0.5+10)
        be=dimCol(1,2)+20+(hfooting-rec);
    end
    
    if (le-dimCol(1,1))*0.5<((hfooting-rec)*0.5+10)
        le=dimCol(1,1)+20+(hfooting-rec);
    end
    
    dimFoot(j,2)=le;
    dimFoot(j,1)=be;
    
    load_conditions_footings(j,1)=napoyo;
    load_conditions_footings(j,2)=pu*0.001;
    contact_pressures(j)=pu/(be*le);
    load_conditions_footings(j,4)=muy*0.00001;
    load_conditions_footings(j,3)=pu*0.001*emin;
end

% To export design results __________________________________________
nombre_archivo='nrebar_footings_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_footings_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_footings.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_footings.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');
%____________________________________________________________________
optimConvFootings=0;
PlotRebarFootingDesign=0;
areaElemFootingsRebar=zeros(nzapatas,2);
efElemFootingsRebar=zeros(nzapatas,2);

collectionArrangementRebarFootings=[];
collectionDispositionRebarFootings=[];

costElemFootings=zeros(nzapatas,1);
nbarsFootingsCollection=zeros(nzapatas,4);
dimensionFootingCollection=zeros(nzapatas,4);
for i=1:nzapatas
    be=dimFoot(i,1);
    le=dimFoot(i,2);
    
    hfooting=hefootings(i);
    col_apoyo=nodesSupportColumns(2,i);
    
    dimCol=dimensions(col_apoyo,:);
    
    [hmodif,mu_axis,bestDispositionFootings,arrangement_bar_footings,...
    nbars_footings,AcBar,bestCost_elem,Ef_elemBar,list_mr_footings]=...
    isrFootings(puSteelFootings,hfooting,be,le,rec,fc,fy,load_conditions_footings,...
    dimCol,colsSymAsymISR,duct,optimConvFootings,PlotRebarFootingDesign);

    hefootings(i)=hmodif;
    hfooting=hmodif;
    collectionArrangementRebarFootings=[collectionArrangementRebarFootings;
                                        arrangement_bar_footings];
    
    collectionDispositionRebarFootings=[collectionDispositionRebarFootings;
                                        bestDispositionFootings];
                                    
              
    
    nbarsFootingsCollection(i,:)=nbars_footings;
    
    dimensionFootingCollection(i,:)=[be le hfooting rec];
    
    % To wirte .txt file for the exportation of results_______________
    if colsSymAsymISR~="ISR"
        fprintf(fileid_03,'%d\n',arrangement_bar_footings);
    
                              
        for j=1:length(bestDispositionFootings(:,1))
            fprintf(fileid_04,'%.2f %.2f\n',bestDispositionFootings(j,:));
        end          
    
    
        fprintf(fileid_01,'%d %d %d %d\n',nbarsFootingsCollection(i,:));
    
    
        fprintf(fileid_02,'%.2f %.2f %.2f %.2f\n',dimensionFootingCollection(i,:));
    
    end
    %__________________________________________________________________
    areaElemFootingsRebar(i,:)=AcBar;
    efElemFootingsRebar(i,:)=Ef_elemBar;
    costElemFootings(i)=bestCost_elem;
end
fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

areaSteelbeams=zeros(nbeams,1);
for i=1:nbeams
    areaSteelbeams(i)=sum(areaSteelbeams3sec(i,:));
end

%%%%%%%------ WEIGHT OF STRUCTURE ______________________________________
[wsteelColsTotal,pacColsElem,wsteelConcBeamsElem,wsteelConcColsElem,...
wsteelConcFootingsElem,vol_beams,vol_cols,vol_footings,wsteelStructure,...
wconc_structure,wbeams,wcols,wfootings,totalWeightStruc]=WeightStruc...
(elem_cols,elem_beams,lenElem,areaElem,areaSteelbeams,areaElemFootingsRebar,...
hefootings,nbeams,ncols,steelAreacols,nzapatas,dimFoot);

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

% COST OF STRUCTURE ____________________________________________________
totalCostStruc=CostStruc(costSteelbeams,cost_steel_cols,...
costElemFootings,fcbeams,fccols,fcfootings,vol_beams,vol_cols,vol_footings);
