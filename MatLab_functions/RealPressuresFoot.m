function [qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_conditions_col,be,le)

%------------------------------------------------------------------------
% Syntax:
% [qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_conditions_col,be,le)

%------------------------------------------------------------------------
% PURPOSE: To compute the distribution of bending moments over the 
% transversal cross-sections of the footing based on the actions applied 
% through the supporting column.
% 
% OUTPUT: qu01,qu02,qu03,qu04: are the distributed pressure at the
%                              upper-right,upper-left,lower-left,lower-right
%                              corners of the isolated footing in the 
%                              reference system orientation (see documentation)
%
%         qprom:               is the average distributed pressure 
%                              considering the four distributed pressures 
%                              at the corners of the element. This pressure
%                              is to be compared with the maximum 
%                              withstanding contact pressure restriction
%
%
% INPUT:  load_conditions_col: is the vector containing the critical load 
%                              condition that the supporting column imposed 
%                              over the isolated footing (Ton,m)
%
%         be,le:               are the transversal dimensions of the isolated
%                              footing on plan view (cm)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

mux=load_conditions_col(1,3);
muy=load_conditions_col(1,4);
pu=load_conditions_col(1,2);
pu=pu*1000; % kg
mumax_x=mux*100000; % kg-cm
mumax_y=muy*100000; % kg-cm

% action on the X-direction - along be
qu01=pu/(be*le)+6*mumax_x/(be*le^2)+mumax_y/(be^2*le);
qu02=pu/(be*le)+6*mumax_x/(be*le^2)-mumax_y/(be^2*le);

qu03=pu/(be*le)-6*mumax_x/(be*le^2)+mumax_y/(be^2*le);
qu04=pu/(be*le)-6*mumax_x/(be*le^2)-mumax_y/(be^2*le);

qprom=(qu01+qu02+qu03+qu04)/4;

