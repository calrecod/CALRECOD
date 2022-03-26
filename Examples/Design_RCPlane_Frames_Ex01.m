clc 
clear all

nnodes=8;
nbars=8;

%%%%%%%%%%%---------materials features---------%%%%%%%%%

% f'c for each element:
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element's material (function of f'c)
Eelem=zeros(nbars,1);
for i=1:nbars
    Eelem(i)=14000*(fpc(i))^0.5;
end

%%%%%%%%%%---------geometrical features---------%%%%%%%%

dimensions=[55 55; % cross-section dimensions of each element
            55 55;
            25 50;
            55 55;
            60 60;
            30 60;
            25 50;
            40 40];

areaElem=zeros(nbars,1);
inertiaElem=zeros(nbars,1);
for i=1:nbars
    areaElem(i)=dimensions(i,1)*dimensions(i,2);

    inertiaElem(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

%coordinates of each node for each bar
coordxy=[0 -100;
         0 400;
         0 800;
         600 800;
         600 400;
         600 -100;
         1200 400;
         1200 -100]; 
     
% coordinates of the beams end' centroids______________________
coordEndBeams=[0 800 0;
               0 400 0;
               600 400 0];
                 
% coordinates of the columns base' centroids______________________
coord_desplante_cols=[0 -100 0;
                      0 400 0;
                      600 400 0;
                      600 -100 0;
                      1200 -100 0];
                  
% coordinates of the footing base' centroids______________________
%_________________________________________________________________

coord_desplante_zapatas=[0 -100 0;
                        600 -100 0;
                        1200 -100 0];

%___________________________________________________________________

%%%---- Initial-final node of each bar --------------------------%%%

ni=[1;2;3;4;5;2;5;7];
nf=[2;3;4;5;6;5;7;8];

lenElem=zeros(nbars,1);
for i=1:nbars
    lenElem(i)=((coordxy(nf(i),1)-coordxy(ni(i),1))^2+(coordxy(nf(i),2)...
        -coordxy(ni(i),2))^2)^0.5;
end    

%  gdl condicion (gdl restringidos) o pre-escritos
bc=[1 0;
    2 0;
    3 0;
    16 0;
    17 0;
    18 0;
    22 0;
    23 0;
    24 0];

% Topology data
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

fglobal=zeros(3*nnodes,1);

type_elem=[1 "Col"; % ID vector to identify beam and column elements
           2 "Col";
           3 "Beam";
           4 "Col";
           5 "Col";
           6 "Beam";
           7 "Beam";
           8 "Col"];
 
beams_LL=[1 100; % Live Load on beams
          2 100;
          3 100];
      
% To detect which and how many beam and column elemets there are_______
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

 np=7; % number of points of analysis for the computation of mechanic 
       % elements for each structural element
 
pu_beams=38.85; % unit construction cost of reinforcement assembly
cols_sym_asym_isr="Asymmetric";
if cols_sym_asym_isr=="Symmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93]; % symmetric rebar
elseif cols_sym_asym_isr=="Asymmetric"
    pu_cols=[32.23, 32.10, 31.96, 31.96, 31.96, 31.96, 31.96;
          36.0, 36.0, 36.0, 36.0, 36.0, 36.0, 36.0]; % asymmetric rebar

elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[28.93];
end

 pu_steel_footings=26.75;
 
fcbeams=300;
fccols=300;
fc_footing=300;

qadm=2.5; % Admissible bearing load of soil (Kg/cm2)
FS=2.0; % Design Safety Factor for isolated footings
nodos_apoyo_columna=[1 6 8; %apoyo
                    1 5 8]; %elemento columna
nfootings=length(nodos_apoyo_columna(1,:));   
nfloors=2;
floor_elem=[1 1 5 6 7 8; % elements that belong to each floor
            2 2 3 4 0 0]; % (see documentation)


%%%%%------------------ Optimization steel of frames...................

%--------------dinamic análysis (inertial forces).......................

support=[1 "Fixed" "Fixed"; % Type of support at the ends of 
       2 "Fixed" "Fixed";   % each element
       3 "Fixed" "Fixed";
       4 "Fixed" "Fixed";
       5 "Fixed" "Fixed";
       6 "Fixed" "Fixed";
       7 "Fixed" "Fixed";
       8 "Fixed" "Fixed"];
           
% Indicate cracking condition in column sections "Cracked" or "Non-cracked"
condition_cracking="Cracked";
%__________________________________________________________________________

unit_weight_elem=zeros(nbars,2);
for i=1:nbars
    unit_weight_elem(i,2)=0.0024; % kg/cm3
end

vertical_loads=zeros(nbars,2);
for i=1:nbeams
    vertical_loads(elem_beams(i),2)=1.1*beams_LL(i,2);
end

fprintf('\nRC FRAME OPTIMIZATION DESIGN\n\n');

% To consider the self weight load as a distributed load in beams
qbarra_y=zeros(nbars,2);
for i=1:nbeams
    qbarra_y(elem_beams(i),2)=1.1*areaElem(elem_beams(i))*...
        unit_weight_elem(elem_beams(i),2)+1.1*(beams_LL(i,2));
    
end

% To consider self weight of columns as punctual vertical loads on its
% supporting nodes
dof_weight_elem_cols=[2 5 14 17 23; % dof
                      1 2 4  5  8]; % cols

for i=1:length(dof_weight_elem_cols(1,:))
    fglobal(dof_weight_elem_cols(1,i))=-1.1*areaElem(dof_weight_elem_cols(2,i))*...
        unit_weight_elem(dof_weight_elem_cols(2,i),2)*lenElem(dof_weight_elem_cols(2,i));

end

%dof_seismic_forces=[4 7]; % in case there are equivalent sesimic loads 
                          % applied
dof_seismic_forces=[];    % in case there are not seismic lateral loads
                          % applied

%%% To compute lateral equvalent seismic forces with the inverted pendulum
% method ______________________________________________________________
if isempty(dof_seismic_forces)==0
    modes=[1];
    sa=200; % ground pseudo-acceleration
    [fmax_floor,modals]=SeismicLoadsNFloor2DFrame(support,unit_weight_elem,...
                    type_elem,lenElem,areaElem,Eelem,inertiaElem,vertical_loads,...
                    nfloors,floor_elem,modes,sa);
    fglobal(dof_seismic_forces)=fmax_floor;
end

%%%----------------------------------------------------------------------
ductility=3; % required ductility on the element's cross-sections

recxy_cols=[4 4];
% To export results
directionData='C:\Users\luizv\OneDrive\Maestría\CAL-RECOD Documentation\Functions\CALRECOD\DesignResults\';
plotAnalysisResults=1; % this parameter is used to indicate if the 
                       % analyisis' results plots will be required or not
                       % (1) yes, (otherwise) no. (see Documentation)

%%% Structural Optimization Design Processan............................
[total_weight_structure,wsteel_cols_total,pac_prom_cols,sectionRestrictions,...
Mp,final_dimensions,displacements,unit_weight_elem,wsteel_conc_beams_elem,...
wsteel_conc_cols_elem,wsteel_conc_footings_elem,he_footings,dim_zap,...
total_cost_structure,inertiaElem,wsteel_structure]=DesignRCPlaneFrameBCIF...
(nnodes,bc,ni,nf,Edof,coordxy,pu_beams,type_elem,pu_cols,nbars,np,...
lenElem,coord_desplante_cols,fpc,inertiaElem,qadm,FS,nodos_apoyo_columna,...
pu_steel_footings,dimensions,Eelem,fcbeams,fccols,fc_footing,fglobal,...
qbarra_y,areaElem,dof_seismic_forces,floor_elem,cols_sym_asym_isr,...
condition_cracking,ductility,elem_cols,elem_beams,recxy_cols,directionData,...
plotAnalysisResults);

% Additional exportation of results for the use of VISUAL-CALRECOD
nombre_archivo='coord_comienzo_beams.txt';
fileid=fopen([directionData,nombre_archivo],'w+t');
for i=1:nbeams
    fprintf(fileid,'%.2f %.2f %.2f\n',coordEndBeams(i,:));
end
fclose(fileid);

nombre_archivo='coord_desplante_columns.txt';

fileid=fopen([directionData,nombre_archivo],'w+t');
for j=1:ncols
    fprintf(fileid,'%.2f %.2f %.2f\n',coord_desplante_cols(j,:));
end
fclose(fileid);

nombre_archivo='coord_desplante_footings.txt';

fileid=fopen([directionData,nombre_archivo],'w+t');
for i=1:nfootings
    fprintf(fileid,'%.2f %.2f %.2f\n',coord_desplante_zapatas(i,:));
end
fclose(fileid);
