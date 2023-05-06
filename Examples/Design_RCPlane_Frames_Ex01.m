% Design_RCPlane_Frames_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally the reinforcement of a structural plane concrete 
%    frame composed of rectangular beams, rectangular columns and isolated
%    footings, subject to gravity and lateral basal loads in both direc-
%    tions (left and right).
%
%    Note: function DesignRCPlaneFrameBCI is the one used to design
%          optimally the reinforcement on all elements. In case it is 
%          necessary some height dimension of the elements might be 
%          increased
%
%          At the end, if required, the design results data can be
%          exported as .csv files for a more proper assessment through the
%          visualization of the optimal designs in ANSYS SpaceClaim or
%          Dynamo which are part of the Visual CALRECOD library. 
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-03-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
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
elem_cols=[];
elem_beams=[];

nbeams=0;
ncols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        nbeams=nbeams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        ncols=ncols+1;
        elem_cols=[elem_cols,j];
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
dimensions=[30 30; % cross-section dimensions of each element
            30 30;
            25 50;
            30 45;
            30 45;
            30 60;
            25 50;
            30 30];

areaElem=zeros(nbars,1);
inertiaElem=zeros(nbars,1);
for i=1:nbars
    areaElem(i)=dimensions(i,1)*dimensions(i,2);
    inertiaElem(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

% coordinates of each node for each bar
coordxy=[0 -100;
         0 400;
         0 800;
         600 800;
         600 400;
         600 -100;
         1200 400;
         1200 -100]; 

%% Location of elements in the 3D system of reference of Visual-CALRECOD
%% Plane Y-Z
% coordinates of the columns base' centroids:
coordBaseCols=[0 0 -100;
              0 0 400;
              0 600 400;
              0 600 -100;
              0 1200 -100];

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
for i=1:nbars
    lenElem(i,1)=((coordxy(nf(i),1)-coordxy(ni(i),1))^2+(coordxy(nf(i),2)...
        -coordxy(ni(i),2))^2)^0.5;
end  

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
condition_cracking="Cracked";

ductility=3; % required ductility on the element's cross-sections for the
             % design of reinforcement

recxy_cols=[4 4]; % concrete cover
            
%% Rebar data
% Commerically available rebars
RebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
%% Optimization design process of reinforced concrete elements

fprintf('\nRC FRAME OPTIMIZATION DESIGN\n\n');

%% Loads
beams_LL=[1 75; % Live and Dead Loads on beams as distributed forces
          2 77;
          3 71];
      
% Unit weight of concrete - To consider self weight
wco=0.0024;
unit_weight_elem=zeros(nbars,2);
for i=1:nbars
    unit_weight_elem(i,2)=wco; % kg/cm3
end

% To consider the self weight load as a distributed load in beams
qbary=zeros(nbars,2);
for i=1:nbeams
    qbary(elem_beams(i),2)=1.1*areaElem(elem_beams(i))*...
        unit_weight_elem(elem_beams(i),2)+1.1*(beams_LL(i,2));
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
Vx1=3500;
Vx2=4500;
fglobal(dof_seismic_forces)=[Vx1;Vx2];

eRefN=5;eRefV=3;eRefM=3; % these are the element's number to be
                         % taken for reference when ploting the mechanical
                         % elements for each bar
[displacements_right,r_right,Ex_right,Ey_right,es_bars_normal_right,...
es_bars_shear_right,es_bars_moment_right]=PlaneFrameStaticLinearAnalysis(nnodes,...
nbars,Eelem,areaElem,inertiaElem,bc,fglobal,ni,nf,qbary,Edof,np,coordxy,0,...
eRefN,eRefV,eRefM);
 
% Lateral forces in the left direction (F-Left):
Vx1=-3500;
Vx2=-4500;
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
    nelem=elem_beams(j);
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
    
    shear_beams(j,:)=es_bars_shear_right(:,nelem);
    load_conditions_beams(j,1)=es_bars_normal_right(1,nelem);
    load_conditions_beams(j,2)=es_bars_moment_right(1,nelem);
    load_conditions_beams(j,3)=max(abs(es_bars_moment_right(2:6,nelem)));
    load_conditions_beams(j,4)=es_bars_moment_right(7,nelem);
end

% -----------------------------------------------------------------------
% Design loads on columns from the structural analysis
% -----------------------------------------------------------------------
emin=0.000001; % no load is considered on the out-of-plane direction

for i=1:ncols
    nelem=elem_cols(i);
    
    % Design loads of seismic forces to the right ---------------------
    % Load at one column's end
    load_conditions_columns(4*i-3,1)=4*i-3;
    
    axial=es_bars_normal_right(1,nelem);
    load_conditions_columns(4*i-3,2)=axial;
    
    [Mux1]=es_bars_moment_right(1,elem_cols(i));
    Mminy1=axial*emin;
    
    load_conditions_columns(4*i-3,3)=Mux1; % Mux (in-plane)
    load_conditions_columns(4*i-3,4)=Mminy1; % Muy (out-of-plane)
    
    % Load at the other column's end
    load_conditions_columns(4*i-2,1)=4*i-2;
    
    axial=es_bars_normal_right(np,nelem);
    load_conditions_columns(4*i-2,2)=axial;
    
    [Mux2]=es_bars_moment_right(np,elem_cols(i));
    Mminy2=axial*emin;
    
    load_conditions_columns(4*i-2,3)=Mux2; % Mux (in-plane)
    load_conditions_columns(4*i-2,4)=Mminy2; % Muy (out-of-plane)
    
    % Design loads of seismic forces to the left ----------------------
    % Load at one column's end
    load_conditions_columns(4*i-1,1)=4*i-1;
    
    axial=es_bars_normal_left(1,nelem);
    load_conditions_columns(4*i-1,2)=axial;
    
    [Mux3]=es_bars_moment_left(1,elem_cols(i));
    Mminy3=axial*emin;
    
    load_conditions_columns(4*i-1,3)=Mux3; % Mux (in-plane)
    load_conditions_columns(4*i-1,4)=Mminy3; % Muy (out-of-plane)
    
    % Load at the other column's end
    load_conditions_columns(4*i,1)=4*i;
    
    axial=es_bars_normal_left(np,nelem);
    load_conditions_columns(4*i,2)=axial;
    
    [Mux4]=es_bars_moment_left(np,elem_cols(i));
    Mminy4=axial*emin;
    
    load_conditions_columns(4*i,3)=Mux4; % Mux (in-plane)
    load_conditions_columns(4*i,4)=Mminy4; % Muy (out-of-plane)
end

% ---------------------------------------------------------------------
% Design loads on Footings from the structural analysis
% ---------------------------------------------------------------------

for j=1:length(r_right)
    if abs(r_right(j))<abs(r_left(j))
        r_right(j)=r_left(j);
    end

end
reactions=r_right;

%% Unit construction cost of rebar assembly for beams, columns and footings
pu_beams=38.85; % unit construction cost of reinforcement assembly
pu_steel_footings=26.75;

cols_sym_asym_isr="Symmetric"; % To choose which type of rebar design
                               % is required
if cols_sym_asym_isr=="Symmetric" || cols_sym_asym_isr=="Asymmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; 
    
elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[29.1];
end
 
%% Directory route to save the design results (if required)
% Note: if not required just set: directionData=[];
directionData=[];
            
%% Optimal design with symmetrical reinforcement in columns       
[totalWeightStruc,wsteelColsTotal,pacColsElem,Mp,final_dimensions,...
unit_weight_elem,wsteelConcBeamsElem,wsteelConcColsElem,...
wsteelConcFootingsElem,hefootings,dimFoot,totalCostStruc,inertiaElem,...
wsteelStructure]=DesignRCPlaneFrameBCI(pu_beams,pu_cols,lenElem,fpc,...
inertiaElem,qadm,FS,nodes_support_column,pu_steel_footings,dimensions,...
fcbeams,fccols,fc_footing,areaElem,cols_sym_asym_isr,RebarAvailable,...
condition_cracking,ductility,elem_cols,elem_beams,recxy_cols,...
load_conditions_beams,load_conditions_columns,reactions,shear_beams,...
coordBaseCols,coordEndBeams,coordBaseFooting,directionData);

disp(' ')
disp('Final dimensions of the elements (symmetrical rebar in columns)')
disp(final_dimensions)

disp('Total weight of the steel reinforcement in the columns')
disp('with symmetrical rebar (Kg)')
disp(wsteelColsTotal);

y1=[wsteelColsTotal];

disp('Total cost of the structure with symmetrical')
disp('reinforcement (MXN)')
disp(totalCostStruc);

y2=[totalCostStruc];

cols_sym_asym_isr="Asymmetric";
if cols_sym_asym_isr=="Symmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; 
elseif cols_sym_asym_isr=="Asymmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93];
elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[28.93];
end

%% Optimal design with asymmetrical reinforcement in columns
directionData=[];
[totalWeightStruc,wsteelColsTotal,pacColsElem,...
Mp,final_dimensions,unit_weight_elem,wsteelConcBeamsElem,...
wsteelConcColsElem,wsteelConcFootingsElem,hefootings,dimFoot,...
totalCostStruc,inertiaElem,wsteelStructure]=DesignRCPlaneFrameBCI...
(pu_beams,pu_cols,lenElem,fpc,inertiaElem,qadm,FS,nodes_support_column,...
pu_steel_footings,dimensions,fcbeams,fccols,fc_footing,areaElem,...
cols_sym_asym_isr,RebarAvailable,condition_cracking,ductility,...
elem_cols,elem_beams,recxy_cols,load_conditions_beams,...
load_conditions_columns,reactions,shear_beams,coordBaseCols,coordEndBeams,coordBaseFooting,directionData);

disp('Final dimensions of the elements (asymmetrical rebar in columns)')
disp(final_dimensions)

disp('Total weight of the steel reinforcement in the structure')
disp('with asymmetrical reinforcement in columns (Kg)')
disp(wsteelColsTotal);

y1=[y1,wsteelColsTotal];

% COST OF STRUCTURE ------------------------------------------------
% ------------------------------------------------------------------

disp('Total cost of the structure with asymmetrical')
disp('reinforcement (MXN)')
disp(totalCostStruc);

y2=[y2,totalCostStruc];

%% Bar diagrams to compare results
x=categorical({'Symmetrical','Asymmetrical'});
figure(10)
bar(x,y1,'black')
hold on
title('Weight of the reinforcing steel in the columns')
xlabel('Type of reinforcement')
ylabel('Weight of the structure')

figure(11)
bar(x,y2,'blue')
hold on
title('Construction cost of the structure: Concrete+Steel')
xlabel('Type of reinforcement')
ylabel('Cost')
