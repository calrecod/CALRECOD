function plotdiagramCircISR(diagramInteraction,conditions)

%------------------------------------------------------------------------
% Syntax:
% plotdiagramCircISR(diagramInteraction,conditions)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To plot the interaction diagram of a circular reinforced 
% column cross-section.
% 
% INPUT:  conditions:           load conditions in format: nload_conditions
%                               rows and four columns as [nload,Pu,Mu]
%
%         diagramInteraction:   interaction diagram of a circular
%                               columns cross-section computed by using
%                               the function: EffCircCols
%                               (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

npdiag=length(diagramInteraction(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%% Not reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagramInteraction(:,2);
y=diagramInteraction(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagramInteraction(:,4);
yfr=diagramInteraction(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=conditions(:,3);
ycondicion=conditions(:,2);
figure(2)
plot(x,y,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([0 max(abs(x))*1.5 y(1,1) y(npdiag,1)]);
grid on
xlabel('Resistant bending moment')
ylabel('Resistant Axial force')
title('ISR Interaction diagram of a circular column')
grid on

