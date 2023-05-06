function totalCostStruc=CostStruc(costSteelBeams,costSteelCols,...
    costSteelFootings,fcbeams,fccols,fcfootings,volbeams,...
    volcols,volfootings,cost_concrete)

%------------------------------------------------------------------------
% Syntax:
% totalCostStruc=CostStruc(costSteelBeams,costSteelCols,...
%   costSteelFootings,fcbeams,fccols,fcfootings,volbeams,...
%   volcols,volfootings,cost_concrete)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the total construction cost of a reinforced concrete 
% plane frame, considering only reinforcing steel and concrete volumes.
% 
% OUTPUT: totalCostStruc:           is the total construction cost estimated 
%                                   for the structural frame, considering
%                                   steel reinforcement and concrete volumes
%
% INPUT:  totalCostStruc:           Is the total construction cost of the
%                                   structural frame, considering both 
%                                   reinforcing steel and concrete volumes
%
%         costSteelBeams:           is the vector containing the construction
%                                   assembly cost of steel reinforcement
%                                   of each beam
%
%         costSteelCols:            is the vector containing the construction
%                                   assembly cost of steel reinforcement
%                                   of each column
%
%         costSteelFootings:        is the vector containing the construction
%                                   assembly cost of steel reinforcement of
%                                   each isolated footing
%
%         fcbeams,fccols,
%         fcfootings:               are the f'c used for beams, 
%                                   columns and isolated footings
%
%         volbeams,volcols,
%         volfootings:              are the total concrete volumes used for 
%                                   each type of element, by using the 
%                                   function WeightStruc
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

total_cost_steel_structure=sum(costSteelBeams)+sum(costSteelCols)+...
    sum(costSteelFootings);

nconcrete=length(cost_concrete(:,1));
for i=1:nconcrete-1
    if cost_concrete(i,1)==fcbeams
        unit_cost_conc_beams=cost_concrete(i,2);
        break;
    elseif cost_concrete(i,1)<fcbeams && cost_concrete(i+1,1)>fcbeams
        unit_cost_conc_beams=0.5*(cost_concrete(i,2)+cost_concrete(i+1,2));
        break;
    end
end

for i=1:nconcrete-1
    if cost_concrete(i,1)==fccols
        unit_cost_conc_cols=cost_concrete(i,2);
        break;
    elseif cost_concrete(i,1)<fccols && cost_concrete(i+1,1)>fccols
        unit_cost_conc_cols=0.5*(cost_concrete(i,2)+cost_concrete(i+1,2));
        break;
    end
end

for i=1:nconcrete-1
    if cost_concrete(i,1)==fcfootings
        unit_cost_conc_footings=cost_concrete(i,2);
        break;
    elseif cost_concrete(i,1)<fcfootings && cost_concrete(i+1,1)>fcfootings
        unit_cost_conc_footings=0.5*(cost_concrete(i,2)+cost_concrete(i+1,2));
        break;
    end
end
cost_conc_beams=volbeams*unit_cost_conc_beams;
cost_conc_cols=volcols*unit_cost_conc_cols;
cost_conc_footings=volfootings*unit_cost_conc_footings;

total_cost_concrete_structure=cost_conc_beams+...
                              cost_conc_cols+...
                              cost_conc_footings;

totalCostStruc=total_cost_concrete_structure+total_cost_steel_structure;
