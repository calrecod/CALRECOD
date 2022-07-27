% ProblemCasesFrames_Ex01_Weight
%----------------------------------------------------------------
% PURPOSE 
%    To determine the total weight of a structural reinforced concrete 
%    plane frame considering concrete volumes and reinforcing steel
%    volumes. The structural frame is composed of rectangular beams,
%    rectangular columns and rectangular/square isolated footings. 
%
%    Note: an avarage of steel percentage area quantity (pmix+pmax)/2 is
%          considered for all element's cross-sections
%
%          function WeightStruc is the one used such weight of the
%          structure
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc 
clear all

nnodes=8;
nbars=8;

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
     
%%%---- Initial-final node of each bar --------------------------%%%

ni=[1;2;3;4;5;2;5;7];
nf=[2;3;4;5;6;5;7;8];

lenElem=zeros(nbars,1);
for i=1:nbars
    lenElem(i)=((coordxy(nf(i),1)-coordxy(ni(i),1))^2+(coordxy(nf(i),2)...
        -coordxy(ni(i),2))^2)^0.5;
end    


type_elem=[1 "Col"; % ID vector to identify beam and column elements
           2 "Col";
           3 "Beam";
           4 "Col";
           5 "Col";
           6 "Beam";
           7 "Beam";
           8 "Col"];
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

hfootings=[30;35;30];
dimFootings=[100 150;
             120 140;
             100 100];

nfootings=length(hfootings);

fcbeams=280;
fdpcbeams=0.85*fcbeams;
fcCols=300;
fcFoot=300;
fdpcFoot=0.85*fcFoot;

fy=4200;

if fcFoot<=280 
    beta1=0.85;
elseif fcFoot>280
    beta1=1.05-fcFoot/1400;
    if beta1<0.65
        beta1=0.65;
    elseif beta1>0.85
        beta1=0.85;
    end
end

pminBeams=(0.7*sqrt(fdpcbeams)/fy);
pmaxBeams=0.025;

pminCols=0.01;
pmaxCols=0.04;

pminFoot=(0.7*sqrt(fdpcFoot)/fy);
pmaxFoot=(0.75*(6000*beta1)/(fy+6000));

areaBarbeams=areaElem(elem_beams)*(pminBeams+pmaxBeams)*0.5;
steelareaCols=areaElem(elem_cols)*(pmaxCols+pminCols)*0.5;

areaBarFootings(:,1)=(dimFootings(:,1).*hfootings+dimFootings(:,2).*hfootings)*...
                     (pminFoot+pmaxFoot)*0.5;
            
areaBarFootings(:,2)=(dimFootings(:,2).*hfootings)*(pminFoot+pmaxFoot)*0.5;
            
[wsteelColsTotal,pacColsElem,wsteelConcBeams,wsteelConcCols,...
    wsteelConcFootings,volbeams,volcols,volfoot,wsteelStruc,wconcStruc,wbeams,...
    wcols,wfootings,weightStructure]=WeightStruc(elem_cols,elem_beams,...
    lenElem,areaElem,areaBarbeams,areaBarFootings,hfootings,nbeams,...
    ncols,steelareaCols,nfootings,dimFootings)

