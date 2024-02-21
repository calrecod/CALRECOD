function UC=unitCostCardColsRec(ucWire,ucRebar,ucWorkers,performWorkers,...
                                unitquanWire,unitquanStirrup,wasteRebar)
%-------------------------------------------------------------------------
% Syntax:
% UC=unitCostCardColsRec(ucWire,ucRebar,ucWorkers,performWorkers,...
%                       unitquanWire,unitquanStirrup,wasteRebar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
%-------------------------------------------------------------------------
% PURPOSE: To compute the unit cost of reinforcment assembly for a
% rectangular concrete column.
% 
% OUTPUT: UC:                   unit cost per unit-weight of reinforcement 
%                               $/Weight
%
% INPUT:  ucWire:               unit cost of wire $/Weight
%
%         ucRebar:              unit cost of rebars $/Weight
%
%         ucWorkers:            unit cost of work labour $/WorkDay
%
%         performWorkers:       performance of labour for reinforcement
%                               assembly Weight/WorkDay
%
%         unitquanWire          Portion of wire per unit-weight of whole 
%                               reinforcement assembly. Common value equal
%                               to 0.04.
%
%         unitquanStirrup       Portion of stirrups per unit-weight of
%                               whole reinforcement assembly. Common value
%                               equal to 0.105.
%
%         wasteRebar            Percentage of waste of rebar per unit- 
%                               weight. A common value is 7%.
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
                                    
sCostWire=ucWire*unitquanWire; % UNITS: $/Kg
sCostStirrups=ucRebar*unitquanStirrup; % UNITS: $/Kg
sCostRebar=ucRebar*(1+wasteRebar*0.01); % UNITS: $/Kg

sManpower=ucWorkers/performWorkers; % UNITS: $/Kg

sTools=sManpower*0.03; % UNITS: $/Kg

UC=sCostWire+...
   sCostStirrups+...
   sCostRebar+...
   sManpower+...
   sTools; % UNITS: $/Kg
