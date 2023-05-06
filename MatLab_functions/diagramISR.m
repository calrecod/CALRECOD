function diagramISR(diagramInteraction,conditions)

%------------------------------------------------------------------------
% Syntax:
% diagramISR(diagramInteraction,conditions)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
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
xcondicion=abs(conditions(:,3));
ycondicion=conditions(:,2);
figure(5)
plot(x,y,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on

plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([0 x(1+fix(npuntos/2)+1)*2 y(1,1) y(npuntos,1)]);
grid on
xlabel('Bending moments')
ylabel('Axial force')
title('Interaction diagram ISR - X')
grid on

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
xcondicion=abs(conditions(:,4));
ycondicion=conditions(:,2);
figure(6)
plot(x,y,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([0 x(1+fix(npuntos/2)+1)*2 y(1,1) y(npuntos,1)]);
grid on
xlabel('Bending moments')
ylabel('Axial force')
title('Interaction diagram ISR - Y')
grid on
