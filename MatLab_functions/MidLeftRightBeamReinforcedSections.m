function MidLeftRightBeamReinforcedSections(h,b,disposition_rebarMid,...
    disposition_rebarLeft,disposition_rebarRight,barTypes1Mid,barTypes2Mid,...
    barTypes1Left,barTypes2Left,barTypes1Right,barTypes2Right)

%------------------------------------------------------------------------
% Syntax:
% MidLeftRightBeamReinforcedSections(h,b,disposition_rebarMid,...
%   disposition_rebarLeft,disposition_rebarRight,barTypes1Mid,barTypes2Mid,...
%   barTypes1Left,barTypes2Left,barTypes1Right,barTypes2Right)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To plot the reinforcement of middle-span, left-span and right-
% span designed beam cross-sections.
% 
% INPUT:  barTypes1Mid,barTypes2Mid:    Vectors that contain the type of 
%                                       rebar of the middle-span cross-
%                                       section of a beam element for both
%                                       in tension and compression, 
%                                       respectively. The vectors size is
%                                       of one column with nrebar rows
%                                       containing a number between 1 and
%                                       7 according to the available
%                                       commercial rebar types stated by
%                                       default
%
%         barTypes1Left,barTypes2Left:  Vectors that contain the type of 
%                                       rebar of the left-span cross-section
%                                       of a beam element for both in tension
%                                       and compression, respectively. The
%                                       vectors size is of one column with
%                                       nrebar rows containing a number 
%                                       between 1 and 7 according to the 
%                                       available commercial rebar types 
%                                       stated by default
%
%         barTypes1Right,
%         barTypes2Right:               Vectors that contain the type of
%                                       rebar of the right-span cross-section
%                                       of a beam element for both in tension
%                                       and compression, respectively. The 
%                                       vectors size is of one column with 
%                                       nrebar rows containing a number
%                                       between 1 and 7 according to the
%                                       available commercial rebar types 
%                                       stated by default
%
%         disposition_rebarMid:         local coordinates of rebars over 
%                                       the middle span cross-section
%
%         disposition_rebarLeft:        local coordinates of rebars over 
%                                       the left span cross-section
%
%         disposition_rebarRight:       local coordinates of rebars over 
%                                       the right span cross-section
%
%         b,h:                          cross-section dimensions
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%% section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordenadas_esq_seccion=[0.5*b 0.5*h;
                         0.5*b -0.5*h;
                         -0.5*b -0.5*h;
                         -0.5*b 0.5*h;
                         0.5*b 0.5*h];
                     
x=coordenadas_esq_seccion(:,1);
y=coordenadas_esq_seccion(:,2);

%% Reinforcement
%-----------------------Central cross-section-------------------------
nbars1=length(barTypes1Mid);
nbars2=length(barTypes2Mid);

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

nbars=length(disposition_rebarMid(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1Mid(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2Mid(i)];
end

for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,disposition_rebarMid(j,1)];
        dispVar1y=[dispVar1y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,disposition_rebarMid(j,1)];
        dispVar2y=[dispVar2y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,disposition_rebarMid(j,1)];
        dispVar3y=[dispVar3y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,disposition_rebarMid(j,1)];
        dispVar4y=[dispVar4y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,disposition_rebarMid(j,1)];
        dispVar5y=[dispVar5y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,disposition_rebarMid(j,1)];
        dispVar6y=[dispVar6y,disposition_rebarMid(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,disposition_rebarMid(j,1)];
        dispVar7y=[dispVar7y,disposition_rebarMid(j,2)];
    end
end
        
figure(2)
plot(x,y,'k')
hold on
xlabel('width')
ylabel('height')
title('Middle Beam Cross-Section')
legend('Cross-section borders')
axis([-(b+20) b+20 -(h+5) h+5])

if isempty(t1)~=1
    figure(2)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor',...
        'red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(2)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor',...
        'blue','DisplayName','Bar Type 5');
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

%-----------------------Left cross-section---------------------------

nbars1=length(barTypes1Left);
nbars2=length(barTypes2Left);

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

nbars=length(disposition_rebarLeft(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1Left(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2Left(i)];
end

for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,disposition_rebarLeft(j,1)];
        dispVar1y=[dispVar1y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,disposition_rebarLeft(j,1)];
        dispVar2y=[dispVar2y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,disposition_rebarLeft(j,1)];
        dispVar3y=[dispVar3y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,disposition_rebarLeft(j,1)];
        dispVar4y=[dispVar4y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,disposition_rebarLeft(j,1)];
        dispVar5y=[dispVar5y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,disposition_rebarLeft(j,1)];
        dispVar6y=[dispVar6y,disposition_rebarLeft(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,disposition_rebarLeft(j,1)];
        dispVar7y=[dispVar7y,disposition_rebarLeft(j,2)];
    end
end
        
figure(1)
plot(x,y,'k')
hold on
xlabel('width')
ylabel('height')
title('Left Beam Cross-Section')
legend('Cross-section borders')
axis([-(b+20) b+20 -(h+5) h+5])

if isempty(t1)~=1
    figure(1)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor',...
        'red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(1)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor',...
        'blue','DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(1)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(1)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(1)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor',...
        'black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(1)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor',...
        'magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(1)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.255 0.069 0]','DisplayName','Bar Type 12');
end

% ----------------------Right cross-section--------------------------

nbars1=length(barTypes1Right);
nbars2=length(barTypes2Right);

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

nbars=length(disposition_rebarRight(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1Right(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2Right(i)];
end

for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,disposition_rebarRight(j,1)];
        dispVar1y=[dispVar1y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,disposition_rebarRight(j,1)];
        dispVar2y=[dispVar2y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,disposition_rebarRight(j,1)];
        dispVar3y=[dispVar3y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,disposition_rebarRight(j,1)];
        dispVar4y=[dispVar4y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,disposition_rebarRight(j,1)];
        dispVar5y=[dispVar5y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,disposition_rebarRight(j,1)];
        dispVar6y=[dispVar6y,disposition_rebarRight(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,disposition_rebarRight(j,1)];
        dispVar7y=[dispVar7y,disposition_rebarRight(j,2)];
    end
end
        
figure(3)
plot(x,y,'k')
hold on
xlabel('width')
ylabel('height')
title('Right Beam Cross-Section')
legend('Cross-section borders')
axis([-(b+20) b+20 -(h+5) h+5])

if isempty(t1)~=1
    figure(3)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor',...
        'red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(3)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor',...
        'blue','DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(3)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(3)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(3)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor',...
        'black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(3)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor',...
        'magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(3)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.255 0.069 0]','DisplayName','Bar Type 12');
end
