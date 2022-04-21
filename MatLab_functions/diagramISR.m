function diagramISR(diagramInteraction,conditions)

%------------------------------------------------------------------------
% Syntax:
% diagramISR(diagramInteraction,conditions)
%
%------------------------------------------------------------------------
% PURPOSE: To graph the interaction diagram of a reinforced column
% cross-section.
% 
% INPUT:  conditions:           load conditions in format: nload_conditions
%                               rows and four columns as [nload,Pu,Mux,Muy]
%
%         diagramInteraction:   interaction diagram data computed by using
%                               the function: widthEfficiencyCols
%                               (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

npuntos=length(diagramInteraction(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% X direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% No reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagramInteraction(:,2);
y=diagramInteraction(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagramInteraction(:,4);
yfr=diagramInteraction(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=conditions(:,3);
ycondicion=conditions(:,2);
figure(5)
plot(x,y,'k')
hold on
plot(xfr,yfr,'r')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red' )
axis([0 x(1+fix(npuntos/2))*2 y(1,1) y(npuntos,1)]);

xlabel('Momentos de flexión')
ylabel('Fuerza axial')
title('Diagrama de Interacción ISR - X')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% No reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagramInteraction(:,6);
y=diagramInteraction(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagramInteraction(:,8);
yfr=diagramInteraction(:,7);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=conditions(:,4);
ycondicion=conditions(:,2);
figure(6)
plot(x,y,'k')
hold on
plot(xfr,yfr,'r')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red' )
axis([0 x(1+fix(npuntos/2))*2 y(1,1) y(npuntos,1)]);

xlabel('Momentos de flexión')
ylabel('Fuerza axial')
title('Diagrama de Interacción ISR - Y')
