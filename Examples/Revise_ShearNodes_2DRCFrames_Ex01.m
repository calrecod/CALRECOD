% Revise_ShearNodes_2DRCFrames_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    1.To design optimally the reinforcement of a structural plane concrete 
%    frame composed of rectangular beams, rectangular columns subject to 
%    gravity and lateral basal loads in both directions (left and right).
%
%    2.To revise their beam-column nodes after the rebar design against 
%    shear forces.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc 
clear all

nnodes=8; % number of nodes of the plane frame
nbars=8; % number of elements of the plane frame

%% Materials and mechanical properties
type_elem=[1 "Col"; % ID vector to identify beam and column elements
           2 "Col";
           3 "Beam";
           4 "Col";
           5 "Col";
           6 "Beam";
           7 "Beam";
           8 "Col"];
       
% f'c for each element type:
fcbeams=250;
fccols=250;
fc_footing=250;

% To detect which how many beam and column elemets there are:
elemcols=[];
elembeams=[];

nbeams=0;
ncols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        nbeams=nbeams+1;
        elembeams=[elembeams,j];
    elseif type_elem(j,2)=="Col"
        ncols=ncols+1;
        elemcols=[elemcols,j];
    end
end

for i=1:nbars
    if type_elem(i,2)=="Col"
        fpc(i,1)=fccols;
    elseif type_elem(i,2)=="Beam"
        fpc(i,1)=fcbeams;
    end
end

% Elasticity modulus of each element's material (function of f'c)
Eelem=zeros(nbars,1);
for i=1:nbars
    Eelem(i)=14000*(fpc(i))^0.5;
end

%% Geometry
dimensions=[30 40; % cross-section dimensions of each element
            30 40;
            25 50;
            30 60;
            30 45;
            30 45;
            25 50;
            30 40];

areaElem=zeros(nbars,1);
inertiaElem=zeros(nbars,1);
for i=1:nbars
    areaElem(i)=dimensions(i,1)*dimensions(i,2);
    inertiaElem(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

% coordinates of each node for each bar
coordxy=[0 0;
         0 400;
         0 800;
         600 800;
         600 400;
         600 0;
         1200 400;
         1200 0]; 

%% Location of elements in the 3D system of reference
%% X-Y-Z
% coordinates of the columns base' centroids:
coordBaseCols=[0 0 0;
              0 0 400;
              0 600 400;
              0 600 0;
              0 1200 0];

% coordinates of the footing base' centroids:

coordBaseFooting=[0 0 -100;
                0 600 -100;
                0 1200 -100];

% coordinates of the beams' starting nodes:
coordEndBeams=[0 0 800;
               0 0 400;
               0 600 400];
           
%% Topology
% Initial-final node of each bar
ni=[1;2;3;4;5;2;5;7];
nf=[2;3;4;5;6;5;7;8];

% Topology conectivity matrix
Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end  

% Length of each element

lenElem(:,1)=((coordxy(nf,1)-coordxy(ni,1)).^2+(coordxy(nf,2)...
        -coordxy(ni,2)).^2).^0.5;

%% Boundary conditions: prescribed dof
bc=[1 0;
    2 0;
    3 0;
    16 0;
    17 0;
    18 0;
    22 0;
    23 0;
    24 0];

%% Additional design and analysis parameters
np=7; % number of points of analysis for the computation of mechanic 
       % elements for each structural element with the FEM

qadm=2.5; % Admissible bearing load of soil (Kg/cm2)
FS=1.5; % Design Safety Factor for isolated footings
nodes_support_column=[1 6 8; % support (dof)
                    1 5 8]; % column element number
                
% Indicate cracking condition in column sections "Cracked" or "Non-cracked"
conditionCracking="Cracked";

duct=3; % required ductility on the element's cross-sections for the
             % design of reinforcement

recCols=[4 4]; % concrete cover
            
%% Rebar data
% Commerically available rebars
rebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
dataCFA=[0,1,1,2]; % Data required for the computation of the 
                   % constructability factor for the rebar designs of columns
                   % Lower and upper values for the range of the CFA and 
                   % weight factors for the uniformity of rebars and rebar 
                   % diameters
                % rebar designs in columns
%% Optimization design process of reinforced concrete elements

fprintf('\nRC FRAME OPTIMIZATION DESIGN\n\n');

%% Loads
beams_LL=[1 30; % Live and Dead Loads on beams as distributed forces
          2 30;
          3 25];
      
% Unit weight of concrete - To consider self weight
wco=0.0024;
unit_weight_elem=zeros(nbars,2);
for i=1:nbars
    unit_weight_elem(i,2)=wco; % kg/cm3
end

% To consider the self weight load as a distributed load in beams
qbary=zeros(nbars,2);
for i=1:nbeams
    qbary(elembeams(i),2)=1.1*areaElem(elembeams(i))*...
        unit_weight_elem(elembeams(i),2)+1.1*(beams_LL(i,2));
end

% To consider self weight of columns as punctual vertical loads on its
% supporting nodes
dofWeightElemCols=[2 5 14 17 23; % dof
                   1 2 4  5  8]; % cols

fglobal=zeros(3*nnodes,1);
for i=1:length(dofWeightElemCols(1,:))
    fglobal(dofWeightElemCols(1,i))=-1.1*areaElem(dofWeightElemCols(2,i))*...
    unit_weight_elem(dofWeightElemCols(2,i),2)*lenElem(dofWeightElemCols(2,i));
end

% To consider lateral forces (seismic or impacto, etc.)
dof_seismic_forces=[4 7]; % in case there are equivalent sesimic loads 
                          % applied

%% Structural Analysis taking all loads in consideration.
% Lateral forces in the right direction (F-Right):                          
Vx1=1500;
Vx2=2000;
fglobal(dof_seismic_forces)=[Vx1;Vx2];

eRefN=5;eRefV=3;eRefM=3; % these are the element's number to be
                         % taken for reference when ploting the mechanical
                         % elements for each bar
[displacements_right,r_right,Ex_right,Ey_right,es_bars_normal_right,...
es_bars_shear_right,es_bars_moment_right]=PlaneFrameStaticLinearAnalysis(nnodes,...
nbars,Eelem,areaElem,inertiaElem,bc,fglobal,ni,nf,qbary,Edof,np,coordxy,0,...
eRefN,eRefV,eRefM);
 
% Lateral forces in the left direction (F-Left):
Vx1=-1500;
Vx2=-2000;
fglobal(dof_seismic_forces)=[Vx1;Vx2];

[displacements_left,r_left,Ex_left,Ey_left,es_bars_normal_left,...
es_bars_shear_left,es_bars_moment_left]=PlaneFrameStaticLinearAnalysis(nnodes,...
nbars,Eelem,areaElem,inertiaElem,bc,fglobal,ni,nf,qbary,Edof,np,coordxy,0,...
eRefN,eRefV,eRefM);

%% Extracting design forces for each element
% -----------------------------------------------------------------------
% Design loads on beams from the structural analysis
% -----------------------------------------------------------------------
for j=1:nbeams
    nelem=elembeams(j);
    for i=1:np
        if abs(es_bars_normal_right(i,nelem))<abs(es_bars_normal_left(i,nelem))
            es_bars_normal_right(i,nelem)=es_bars_normal_left(i,nelem);

        end
        if abs(es_bars_shear_right(i,nelem))<abs(es_bars_shear_left(i,nelem))
            es_bars_shear_right(i,nelem)=es_bars_shear_left(i,nelem);

        end
        if abs(es_bars_moment_right(i,nelem))<abs(es_bars_moment_left(i,nelem))
            es_bars_moment_right(i,nelem)=es_bars_moment_left(i,nelem);

        end
    end
    
    load_conditions_beams(j,1)=es_bars_normal_right(1,nelem);
    load_conditions_beams(j,2)=es_bars_moment_right(1,nelem);
    load_conditions_beams(j,3)=-max(abs(es_bars_moment_right(2:6,nelem)));
    load_conditions_beams(j,4)=es_bars_moment_right(7,nelem);
end

% -----------------------------------------------------------------------
% Design loads on columns from the structural analysis
% -----------------------------------------------------------------------
emin=0.000001; % no load is considered on the out-of-plane direction

for i=1:ncols
    nelem=elemcols(i);
    
    % Design loads of seismic forces to the right ---------------------
    % Load at one column's end
    load_conditions_columns(4*i-3,1)=4*i-3;
    
    axial=es_bars_normal_right(1,nelem);
    load_conditions_columns(4*i-3,2)=axial;
    
    [Mux1]=es_bars_moment_right(1,elemcols(i));
    Mminy1=axial*emin;
    
    load_conditions_columns(4*i-3,3)=Mux1; % Mux (in-plane)
    load_conditions_columns(4*i-3,4)=Mminy1; % Muy (out-of-plane)
    
    % Load at the other column's end
    load_conditions_columns(4*i-2,1)=4*i-2;
    
    axial=es_bars_normal_right(np,nelem);
    load_conditions_columns(4*i-2,2)=axial;
    
    [Mux2]=es_bars_moment_right(np,elemcols(i));
    Mminy2=axial*emin;
    
    load_conditions_columns(4*i-2,3)=Mux2; % Mux (in-plane)
    load_conditions_columns(4*i-2,4)=Mminy2; % Muy (out-of-plane)
    
    % Design loads of seismic forces to the left ----------------------
    % Load at one column's end
    load_conditions_columns(4*i-1,1)=4*i-1;
    
    axial=es_bars_normal_left(1,nelem);
    load_conditions_columns(4*i-1,2)=axial;
    
    [Mux3]=es_bars_moment_left(1,elemcols(i));
    Mminy3=axial*emin;
    
    load_conditions_columns(4*i-1,3)=Mux3; % Mux (in-plane)
    load_conditions_columns(4*i-1,4)=Mminy3; % Muy (out-of-plane)
    
    % Load at the other column's end
    load_conditions_columns(4*i,1)=4*i;
    
    axial=es_bars_normal_left(np,nelem);
    load_conditions_columns(4*i,2)=axial;
    
    [Mux4]=es_bars_moment_left(np,elemcols(i));
    Mminy4=axial*emin;
    
    load_conditions_columns(4*i,3)=Mux4; % Mux (in-plane)
    load_conditions_columns(4*i,4)=Mminy4; % Muy (out-of-plane)
end

%% Unit construction cost of rebar assembly for beams, columns and footings
puBeams=38.85; % unit construction cost of reinforcement assembly
pu_steel_footings=26.75;

colsSymAsymISR="Symmetric"; % To choose which type of rebar design
                               % is required

puColsISR=[29.1];

%% ----------------------------------------------------------------------
%                  OPTIMAL DESIGN OF REBAR ON ELEMENTS
% -----------------------------------------------------------------------

%% DESIGN OF BEAMS 
nbeams=length(elembeams);
                
h_rec_beams=4; % concrete cover for beams
b_rec_beams=4; % (cm)                 
fy=4200; % yield stress of steel reinforcement
wac=7.8e-3; % unit volume weight of the reinforcement (Kg/cm3)
wco=2.4e-3; % unit volume weight of the concrete (Kg/cm3)
graphConvergenceBeamPlot=0;
rebarBeamDesignPlots=0;

for i=1:nbeams
    nelem=elembeams(i);
    b=dimensions(nelem,1); %cm
    h=dimensions(nelem,2); %cm
    span=lenElem(nelem); %cm
    fc=fpc(nelem);

    [sepbarsRestric,b,h,inertia_modif,dispositionBar_Der,barArrangementDerComp,...
    barArrangementDerTens,dispositionBar_Center,barArrangementCentralTens,...
    barArrangementCentralComp,dispositionBar_Izq,barArrangementIzqTens,...
    barArrangementIzqComp,minAreabar3sec,Efelem3secISR,bestCostVar,efrebar3sec,...
    minAreaRebarAverage,Mr3section]=beamsISR(puBeams,span,wac,b,h,h_rec_beams,...
    rebarAvailable,fc,fy,load_conditions_beams(i,:),colsSymAsymISR,duct,b_rec_beams,...
    rebarBeamDesignPlots,graphConvergenceBeamPlot);

    % To collect and save design results:
    MrbeamsCollection(i,:)=Mr3section;
                             
    dim_beams_collection(i,:)=[b h span h_rec_beams];
             
    % Update dimensions of beams in the global model (in case they were
    % modified during the design process):
    dimensions(nelem,2)=h; %cm 
    
    areaSteelbeams3sec(i,:)=minAreabar3sec;

    % To compute the resistance at each end of the beam elements
    Mp(elembeams(i),1)=MrbeamsCollection(i,1);
    Mp(elembeams(i),2)=MrbeamsCollection(i,3);
end

%% DESIGN OF COLUMNS
ncols=length(elemcols);

optimaColConvPlot=0;
plotISRColResults=0;
plotRebarDesign=0;

puCostCardBuild=[1.2,0.9,128,214.75,0.04,0.105,7];
nlf=length(load_conditions_columns(:,1))/ncols; % This variable is to
                                                % to know how many loads
                                                % take place for each 
                                                % column.
for i=1:ncols
    nelem=elemcols(i);
    b=dimensions(nelem,1); %cm
    h=dimensions(nelem,2); %cm
    
    height=lenElem(nelem); %cm
    fc=fpc(nelem); %kg/cm2
    
    loadCol=load_conditions_columns(nlf*i-(nlf-1):nlf*i,:);
    
    [Inertia_xy_modif,b,h,bestArrangement,bestdisposicionRebar,costCol,...
    AcSecCols,EfsecCol,Mr_col,bestCFA]=isrColumnsSymAsym(puColsISR,height,...
    b,h,recCols,fy,fc,loadCol,colsSymAsymISR,conditionCracking,duct,wac,...
    rebarAvailable,puCostCardBuild,dataCFA,optimaColConvPlot,...
    plotISRColResults,plotRebarDesign);
    
    MrColumnsCollection(i,:)=Mr_col;
    
    dimensions(nelem,2)=h; % in case hte height dimension was modified 
    Mp(elemcols(i),1)=MrColumnsCollection(i,1);
    Mp(elemcols(i),2)=MrColumnsCollection(i,1);
end

%% NODE REVISION
% Shear resistance
[VcrxNode,VuxNode,VxEfNodes,hnodeMin]=shearNodes2DFramesSI(ni,nf,...
 elembeams,elemcols,fccols,fy,MrbeamsCollection,areaSteelbeams3sec,...
 dimensions,nnodes,coordBaseCols,lenElem);

disp('Final dimensions of the elements (symmetrical rebar in columns)')
disp(dimensions)
