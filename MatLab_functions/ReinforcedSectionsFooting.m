function ReinforcedSectionsFooting(h,be,le,dispositionRebar1,...
    dispositionRebar2,barTypes1B,barTypes2B,barTypes1L,barTypes2L)

%------------------------------------------------------------------------
% Syntax:
% ReinforcedSectionsFooting(h,be,le,dispositionRebar1,...
%   dispositionRebar2,barTypes1B,barTypes2B,barTypes1L,barTypes2L)
%
%------------------------------------------------------------------------
% PURPOSE: To plot both reinforced concrete sections of a rectangular
% isolated footing.
% 
% INPUT:  h:                    footing height
%
%         be,le:                transversal cross-section dimensions 
%
%         dispositionRebar1,
%         dispositionRebar2:    local rebar coordinates over the transversal
%                               cross-section with the $be$ and le dimension,
%                               respectively
%
%         barTypes1B,
%         barTypes1L:           types of rebar in tension for both 
%                               transversal cross-sections, be and le 
%                               dimension, respectively
%
%         barTypes2B,
%         barTypes2L:           types of rebar in compression for both
%                               transversal cross-sections, be and le 
%                               dimension, respectively
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%% section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_var=dispositionRebar1(:,1);
y_var=dispositionRebar1(:,2);

%%---------------------------------------------%%%
coordenadas_esq_seccion=[0.5*be 0.5*h;
                         0.5*be -0.5*h;
                         -0.5*be -0.5*h;
                         -0.5*be 0.5*h;
                         0.5*be 0.5*h];
                     
x=coordenadas_esq_seccion(:,1);
y=coordenadas_esq_seccion(:,2);

nbars1=length(barTypes1B);
nbars2=length(barTypes2B);

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

nbars=length(dispositionRebar1(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1B(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2B(i)];
end


for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,dispositionRebar1(j,1)];
        dispVar1y=[dispVar1y,dispositionRebar1(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,dispositionRebar1(j,1)];
        dispVar2y=[dispVar2y,dispositionRebar1(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,dispositionRebar1(j,1)];
        dispVar3y=[dispVar3y,dispositionRebar1(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,dispositionRebar1(j,1)];
        dispVar4y=[dispVar4y,dispositionRebar1(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,dispositionRebar1(j,1)];
        dispVar5y=[dispVar5y,dispositionRebar1(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,dispositionRebar1(j,1)];
        dispVar6y=[dispVar6y,dispositionRebar1(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,dispositionRebar1(j,1)];
        dispVar7y=[dispVar7y,dispositionRebar1(j,2)];
    end
end
        
figure(3)
plot(x,y,'k')
hold on
xlabel('Width (cm)')
ylabel('Height (cm)')
title('Footing section - B')
legend('Footing boundaries ')
axis([-(be+20) be+20 -(h+5) h+5])
%gtext('rec=3cm')

if isempty(t1)~=1
    figure(3)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor','red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(3)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor','blue','DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(3)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor','[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(3)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor','[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(3)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor','black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(3)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor','magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(3)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor','[0.255 0.069 0]','DisplayName','Bar Type 12');
end
        
% REINFORCEMENT ON THE OTHER SECTION _________________________________

%%%%%%%%%%%%%%%%%%%%%%%%% seccion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_var=dispositionRebar2(:,1);
y_var=dispositionRebar2(:,2);

%%---------------------------------------------%%%
coordenadas_esq_seccion=[0.5*le 0.5*h;
                         0.5*le -0.5*h;
                         -0.5*le -0.5*h;
                         -0.5*le 0.5*h;
                         0.5*le 0.5*h];
                     
x=coordenadas_esq_seccion(:,1);
y=coordenadas_esq_seccion(:,2);

nbars1=length(barTypes1L);
nbars2=length(barTypes2L);
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

nbars=length(dispositionRebar2(:,1));
barTypes=[];
for i=1:nbars1
    barTypes=[barTypes,barTypes1L(i)];
end

for i=1:nbars2
    barTypes=[barTypes,barTypes2L(i)];
end

for j=1:nbars
    if barTypes(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,dispositionRebar2(j,1)];
        dispVar1y=[dispVar1y,dispositionRebar2(j,2)];
    elseif barTypes(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,dispositionRebar2(j,1)];
        dispVar2y=[dispVar2y,dispositionRebar2(j,2)];
    elseif barTypes(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,dispositionRebar2(j,1)];
        dispVar3y=[dispVar3y,dispositionRebar2(j,2)];
    elseif barTypes(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,dispositionRebar2(j,1)];
        dispVar4y=[dispVar4y,dispositionRebar2(j,2)];
    elseif barTypes(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,dispositionRebar2(j,1)];
        dispVar5y=[dispVar5y,dispositionRebar2(j,2)];
    elseif barTypes(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,dispositionRebar2(j,1)];
        dispVar6y=[dispVar6y,dispositionRebar2(j,2)];
    elseif barTypes(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,dispositionRebar2(j,1)];
        dispVar7y=[dispVar7y,dispositionRebar2(j,2)];
    end
end
        
figure(2)
plot(x,y,'k')
hold on
xlabel('Width (cm)')
ylabel('Height (cm)')
title('Footing section - L')
legend('Footing boundaries ')
axis([-(le+20) le+20 -(h+5) h+5])
%gtext('rec=3cm')

if isempty(t1)~=1
    figure(2)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor','red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(2)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor','blue','DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(2)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor','[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(2)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor','[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(2)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor','black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(2)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor','magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(2)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor','[0.255 0.069 0]','DisplayName','Bar Type 12');
end