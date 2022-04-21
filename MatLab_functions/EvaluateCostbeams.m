function cost = EvaluateCostbeams(nvHor1,nvHor2,arreglo_t1,arreglo_t2,...
                                    pu,availableRebar)

%------------------------------------------------------------------------
% Syntax:
% cost = EvaluateCostbeams(nvHor1,nvHor2,arreglo_t1,arreglo_t2,...
%                                    pu,availableRebar)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the unit linear cost of rebar assembly for a beam 
% cross-section given an average unit cost in units ($/Kg)
% 
% OUTPUT: cost:             unit linear cost of rebar for a beam cross-
%                           section in units ($/cm)
%
% INPUT:  arreglo_t1,       Vectors that contain the type of rebar for the
%         arreglo_t2:       optimal option both in tension and compression,
%                           respectively. The vectors size is of one column 
%                           with nrebar rows containing a number between 1 
%                           and 7 according to the available commercial 
%                           rebar types stated by default
%
%         nvHor1,nvHor2:    number of rebars in tension and compression,
%                           respectively
%
%         pu:               unit-cost of rebar assebmly for a beam
%                           cross-section
%
%         availableRebar:   data base of the available commercial types of 
%                           rebar
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
                        
cost=0;
for i=1:nvHor1
    cost=cost+availableRebar(arreglo_t1(i),5)*pu(1)*0.01;
end

for j=1:nvHor2
    cost=cost+availableRebar(arreglo_t2(j),5)*pu(1)*0.01;
end
 
