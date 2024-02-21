function [qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_col,be,le,...
                                     typeFooting,dimCol,plotPressure)

%------------------------------------------------------------------------
% Syntax:
% [qu01,qu02,qu03,qu04,qprom]=RealPressuresFoot(load_col,be,le,typeFooting,
%                                               dimCol,plotPressure)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
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
%         typeFooting:         1 - Standard Isolated footing
%                              2 - Bordering Isolated footing
%                              3 - Corner Isolated footing
%
%         dimCol:              column cross-section dimensions [b,h]
%
%         plotPressure:        Optional argument to plot the soil's
%                              pressure distribution over the footing:
%                              0 - Do not plot 
%                              1 - Do plot
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

mux=load_col(1,3);
muy=load_col(1,4);
pu=load_col(1,2);

% Computing the column's loads excentricities
if typeFooting==1
    xr=mux/pu+be/2;
    ex=be/2-xr;
    yr=-muy/pu+le/2;   
    ey=yr-le/2;
    
elseif typeFooting==2
    hc=dimCol(2);

    xr=mux/pu+hc/2;
    ex=be/2-xr;  
    
    yr=-muy/pu+le/2;
    ey=yr-le/2;

elseif typeFooting==3
    hc=dimCol(2);
    bc=dimCol(1);
    xr=mux/pu+hc/2;
    ex=be/2-xr;  
    
    yr=-muy/pu+bc/2;
    ey=yr-le/2;
end

% Contact forces on each footing's corner. 
% (+ means that the pressure is in compression with the soil)
qu03=-pu/(be*le)*(1+6*ex/be+6*ey/le);
qu02=-pu/(be*le)*(1+6*ex/be-6*ey/le);

qu04=-pu/(be*le)*(1-6*ex/be+6*ey/le);
qu01=-pu/(be*le)*(1-6*ex/be-6*ey/le);

qui=[qu01,qu02,qu03,qu04];

qprom=0;
for i=1:4
    if qui(i)>=0 % only positive contact pressures are considered
                 % (when the soil is in compression)
        qprom=qprom+qui(i);
    end
end
qprom=qprom/4;
if nargin>5
    if plotPressure==1
        % Visualization of the stresses distribution
        xValues=[-0.5*be 0.5*be;
                -0.5*be 0.5*be];

        yValues=[-0.5*le -0.5*le;
                0.5*le 0.5*le];

        zValues=[qu03 qu04;
                qu02 qu01];

        figure(1)
        surfc(xValues,yValues,zValues);
        colormap default;
        shading interp;
        colorbar
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Distribution of contact pressures over the footing')
    end
end