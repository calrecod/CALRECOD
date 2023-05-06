function beamReinforcedSection(h,b,disposition_rebar,barTypes1,barTypes2)

%------------------------------------------------------------------------
% Syntax:
% beamReinforcedSection(h,b,disposition_rebar,barTypes1,barTypes2)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To plot the reinforcement of a rectangular cross-section.
% 
% INPUT:  barTypes1,barTypes2:          Vectors that contain the type of 
%                                       rebar for the optimal option both
%                                       in tension and compression, 
%                                       respectively. The vectors size is 
%                                       of one column with nrebar rows 
%                                       containing a number between 1 and 7
%                                       according to the available commercial
%                                       rebar types stated by default
%
%         disposition_rebar:            local coordinates of rebars over
%                                       the cross-section
%
%         b,h:                          cross-section dimensions
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%% section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%---------------------------------------------%%%
coordenadas_esq_seccion=[0.5*b 0.5*h;
                         0.5*b -0.5*h;
                         -0.5*b -0.5*h;
                         -0.5*b 0.5*h;
                         0.5*b 0.5*h];
                     
x=coordenadas_esq_seccion(:,1);
y=coordenadas_esq_seccion(:,2);

nbars1=length(barTypes1);
nbars2=length(barTypes2);

%-------------- baem plot ------------------%

t1=[];
t2=[];
t3=[];
t4=[];
t5=[];
t6=[];
t7=[];

dispVar1x=[];
dispVar1y=[];

dispVar2x=[];
dispVar2y=[];

dispVar3x=[];
dispVar3y=[];

dispVar4x=[];
dispVar4y=[];

dispVar5x=[];
dispVar5y=[];

dispVar6x=[];
dispVar6y=[];

dispVar7x=[];
dispVar7y=[];

nbars=length(disposition_rebar(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2(i)];
end

for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,disposition_rebar(j,1)];
        dispVar1y=[dispVar1y,disposition_rebar(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,disposition_rebar(j,1)];
        dispVar2y=[dispVar2y,disposition_rebar(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,disposition_rebar(j,1)];
        dispVar3y=[dispVar3y,disposition_rebar(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,disposition_rebar(j,1)];
        dispVar4y=[dispVar4y,disposition_rebar(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,disposition_rebar(j,1)];
        dispVar5y=[dispVar5y,disposition_rebar(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,disposition_rebar(j,1)];
        dispVar6y=[dispVar6y,disposition_rebar(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,disposition_rebar(j,1)];
        dispVar7y=[dispVar7y,disposition_rebar(j,2)];
    end
end
        
% Central cross-section layout ________________________________________
figure(2)
plot(x,y,'k')
hold on
xlabel('width')
ylabel('height')
title('Rectangular Cross-Section')
legend('Section borders')
axis([-(b+20) b+20 -(h+5) h+5])

if isempty(t1)~=1
    figure(2)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor','red',...
        'DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(2)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor','blue',...
        'DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(2)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(2)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(2)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor',...
        'black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(2)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor',...
        'magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(2)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.255 0.069 0]','DisplayName','Bar Type 12');
end
        
