function diagramISRot(CoordCorners,load_conditions,diagram,b,h)   
%------------------------------------------------------------------------
% Syntax:
% diagramISRot(CoordCorners,load_conditions,diagram,b,h)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To plot a rotated rectangular column cross-section
% as well as its interaction diagram. No rebars are plot given that the
% cross-section is considered to be reinforced with the ISR.
% 
% INPUT:  CoordCorners:         is the vector 4 x 2 containing the four
%                               pair of coordinates of each cross-section 
%                               corner
%
%         load_conditions:      load conditions in format: nload_conditions
%                               rows and four columns as [nload,Pu,Mu]
%
%         diagram:              is the interaction diagram data 
%                               (see Documentation)
%
%         b,h:                  are the cross-section dimensions (width and
%                               heigth, respectively)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

CoordCorners(5,:)=CoordCorners(1,:);                      
x=CoordCorners(:,1);
y=CoordCorners(:,2);

%-------------------------- column plot ------------------------------%

figure(3)
plot(x,y,'k -','linewidth',1)
hold on
xlabel('x´')
ylabel('y´')
title('Rectangular rotated column')
legend('Column´s boundaries')
axis([-(b+10) b+10 -(h+10) h+10])
        
%------------------------ end column plot ----------------------------%

%%%%%%%%%%%%%%%%%%%%%% Interaction diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Not reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagram(:,2);
y=diagram(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagram(:,4);
yfr=diagram(:,3);
npuntos=length(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=abs(load_conditions(:,3));
ycondicion=load_conditions(:,2);
figure(2)
plot(x,y,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.05,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([0 diagram(fix(npuntos/2)+1,2)*2 diagram(1,1) diagram(npuntos,1)]);
xlabel('Flexure moment ')
ylabel('Axial Force')
title(strcat('Interaction diagram - Rotated rectangular column',...
    'reinforced with the ISR'))
grid on