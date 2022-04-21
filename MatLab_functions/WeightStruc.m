function [wsteelColsTotal,pacColsElem,wsteelConcBeams,wsteelConcCols,...
    wsteelConcFootings,volbeams,volcols,volfoot,wsteelStruc,wconcStruc,wbeams,...
    wcols,wfootings,weightStructure]=WeightStruc(elem_cols,elem_beams,...
    lenElem,areaElem,areaBarbeams,areaBarFootings,hfootings,nbeams,...
    ncols,steelareaCols,nfootings,dimFootings)

%------------------------------------------------------------------------
% Syntax:
% [wsteelColsTotal,pacColsElem,wsteelConcBeams,wsteelConcCols,...
%   wsteelConcFootings,volbeams,volcols,volfoot,wsteelStruc,wconcStruc,wbeams,...
%   wcols,wfootings,weightStructure]=WeightStruc(elem_cols,elem_beams,...
%   lenElem,areaElem,areaBarbeams,areaBarFootings,hfootings,nbeams,...
%   ncols,steelareaCols,nfootings,dimFootings)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the weight of a reinforced concrete plane frame 
% structure as well as all of its elements individually. All units are in 
% Kg,cm
% 
% OUTPUT: wsteelColsTotal:          Sum of the total weight of reinforcing
%                                   steel of each column
%
%         pacColsElem:              Vector containing the percentage area
%                                   of steel reinforcement on each of the 
%                                   columns. Size = [ncolumns,1]
%
%         wsteelConcBeams:          Sum of the weight of both steel
%                                   reinforcement and concrete of each of
%                                   the beams
%
%         wsteelConcCols:           Sum of the weight of both steel
%                                   reinforcement and concrete of each column
%
%         wsteelConcFootings:       Sum of the weight of both steel
%                                   reinforcement and concrete of each 
%                                   isolated footing
%
%         volbeams,volcols,volfoot: are sum of the volume of concrete of
%                                   each beam, column and isolated footing, 
%                                   respectively
%
%         wsteelStruc:              Total weight of steel reinforcement of
%                                   the whole structural frame
%
%         wconcStruc:               Total weight of concrete of the whole
%                                   structural frame
%
%         wbeams,wcols,wfootings:   sum of the total weight of each type of
%                                   element; beams, columns and isolated 
%                                   footings, respectively 
%
%         weightStructure:          Total weight of the structure, 
%                                   considering both steel reinforcement 
%                                   and concrete volumes
%
% INPUT:  elem_cols:                is the vector containing the element
%                                   number code of those elements identified
%                                   as columns
%
%         elem_beams:               is the vector containing the element 
%                                   number code of those elements identified
%                                   as beams
%
%         lenElem:                  is a vector containing the length of
%                                   each element
%
%         areaElem:                 is the vector containing the cross-
%                                   section area of each element
%
%         areaBarbeams:             is the vector containing the quantity
%                                   of steel rebar area of each beam element
%
%         areaBarFootings:          is the vector containing the quantity 
%                                   of steel rebar area of each isolated 
%                                   footing element
%
%         hfootings:                is the vector containing the design
%                                   height dimension of each isolated footing
%
%         nbeams,ncols:             are the total number of beam and column
%                                   elements, respectively
%
%         steelareaCols:            is the vector containing the total rebar 
%                                   area for each column element
%
%         nfootings:                is the number of isolate footing elements
%
%         dimFootings:              is the vector containing the transversal
%                                   design cross-sections of each isolated 
%                                   footing
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%% WEIGHT OF STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%

wsteelConcBeams=zeros(nbeams,1);
wconcStruc=0;
wsteelStruc=0;
volbeams=0;
wbeams=0;
for i=1:nbeams
    nelem=elem_beams(i);
    wconc=lenElem(nelem)*areaElem(nelem)*0.000001*2400; %kg
    volbeams=volbeams+lenElem(nelem)*areaElem(nelem)*0.000001; % cm3
    wsteel=areaBarbeams(i)*0.000001*lenElem(nelem)*7800; %kg
    
    wbeams=wbeams+wconc+wsteel;
    wconcStruc=wconcStruc+wconc;
    wsteelStruc=wsteelStruc+wsteel;
    
    wsteelConcBeams(i)=wsteel+wconc;
end

wsteelConcCols=zeros(ncols,1);
wsteelColsTotal=0;
pacColsElem=zeros(ncols,1);
wcols=0;
volcols=0;
for i=1:ncols
    nelem=elem_cols(i);
    wconc=lenElem(nelem)*areaElem(nelem)*0.000001*2400; %kg
    volcols=volcols+lenElem(nelem)*areaElem(nelem)*0.000001;
    
    wsteel=steelareaCols(i)*0.000001*lenElem(nelem)*7800; %kg
    wsteelColsTotal=wsteelColsTotal+wsteel;
    pacColsElem(i)=steelareaCols(i)/(areaElem(nelem));
    
    wsteelConcCols(i)=wsteel+wconc;
    wcols=wcols+wconc+wsteel;
    wconcStruc=wconcStruc+wconc;
    wsteelStruc=wsteelStruc+wsteel;
end
pac_prom_cols=sum(pacColsElem)/ncols;

wsteelConcFootings=zeros(nfootings,1);
wfootings=0;
volfoot=0;
for i=1:nfootings
    be=dimFootings(i,1);
    le=dimFootings(i,2);
    wconc=be*le*hfootings(i)*0.000001*2400; %kg
    volfoot=volfoot+be*le*hfootings(i)*0.000001;
    
    wsteel=(areaBarFootings(i,1)*0.0001*be*0.01+...
            areaBarFootings(i,2)*0.0001*le*0.01)*7800; %kg
    
    wsteelConcFootings(i)=wsteel+wconc;
    
    wfootings=wfootings+wconc+wsteel;
    wconcStruc=wconcStruc+wconc;
    wsteelStruc=wsteelStruc+wsteel;
end

weightStructure=wfootings+wcols+wbeams;

WeightPercentageConcStruc=wconcStruc/weightStructure;
WeightPercentageSteelStruc=wsteelStruc/weightStructure;

 %%%%%%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%
