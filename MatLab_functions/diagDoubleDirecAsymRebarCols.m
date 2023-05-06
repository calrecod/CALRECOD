function diagDoubleDirecAsymRebarCols(load_conditions,diagrama1,diagrama2,...
            dispositionRebar,h,b,arregloVar)

%------------------------------------------------------------------------
% Syntax:
% diagDoubleDirecAsymRebarCols(load_conditions,diagrama1,diagrama2,
% dispositionRebar,h,b,arregloVar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To graph the main four interaction diagrams of a reinforced 
% column's cross-section (2 interaction diagrams for each cross-section's 
% axis) and the reinforced cross-section layout itself. 
% 
% INPUT:  load_conditions:      load conditions in format: nload_conditions
%                               rows and four columns as [nload,Pu,Mux,Muy]
%
%         diagrama1:            interaction diagram data with positive
%                               bending moments
%
%         diagrama2:            interaction diagram data with negative
%                               bending moments
%
%         dispositionRebar:     rebar local coordinates over the
%                               cross-section [x,y]
%
%         arregloVar:           list of types of rebar distributed over the
%                               cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-08-02
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%% Cross-section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordenadas_esq_seccion=[0.5*b 0.5*h;
                         0.5*b -0.5*h;
                         -0.5*b -0.5*h;
                         -0.5*b 0.5*h;
                         0.5*b 0.5*h];
                     
x=coordenadas_esq_seccion(:,1);
y=coordenadas_esq_seccion(:,2);

nbars=length(arregloVar);

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

for j=1:nbars
    if arregloVar(j)==1
        t1=[t1,1];
        dispVar1x=[dispVar1x,dispositionRebar(j,1)];
        dispVar1y=[dispVar1y,dispositionRebar(j,2)];
    elseif arregloVar(j)==2
        t2=[t2,2];
        dispVar2x=[dispVar2x,dispositionRebar(j,1)];
        dispVar2y=[dispVar2y,dispositionRebar(j,2)];
    elseif arregloVar(j)==3
        t3=[t3,3];
        dispVar3x=[dispVar3x,dispositionRebar(j,1)];
        dispVar3y=[dispVar3y,dispositionRebar(j,2)];
    elseif arregloVar(j)==4
        t4=[t4,4];
        dispVar4x=[dispVar4x,dispositionRebar(j,1)];
        dispVar4y=[dispVar4y,dispositionRebar(j,2)];
    elseif arregloVar(j)==5
        t5=[t5,5];
        dispVar5x=[dispVar5x,dispositionRebar(j,1)];
        dispVar5y=[dispVar5y,dispositionRebar(j,2)];
    elseif arregloVar(j)==6
        t6=[t6,6];
        dispVar6x=[dispVar6x,dispositionRebar(j,1)];
        dispVar6y=[dispVar6y,dispositionRebar(j,2)];
    elseif arregloVar(j)==7
        t7=[t7,7];
        dispVar7x=[dispVar7x,dispositionRebar(j,1)];
        dispVar7y=[dispVar7y,dispositionRebar(j,2)];
    end
end

figure(4)
plot(x,y,'k -','linewidth',1)
hold on
xlabel('Width')
ylabel('Height')
title('Rectangular Column')
legend('Column´s boundaries')
axis([-(b+20) b+20 -(h+5) h+5])

if isempty(t1)~=1
    figure(4)
    plot(dispVar1x,dispVar1y,'r o','linewidth',1,'MarkerFaceColor',...
        'red','DisplayName','Bar Type 4');
end
if isempty(t2)~=1
    figure(4)
    plot(dispVar2x,dispVar2y,'b o','linewidth',1,'MarkerFaceColor',...
        'blue','DisplayName','Bar Type 5');
end
if isempty(t3)~=1
    figure(4)
    plot(dispVar3x,dispVar3y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.05 0.205 0.05]','DisplayName','Bar Type 6');
end
if isempty(t4)~=1
    figure(4)
    plot(dispVar4x,dispVar4y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.072 0.061 0.139]','DisplayName','Bar Type 8');
end
if isempty(t5)~=1
    figure(4)
    plot(dispVar5x,dispVar5y,'k o','linewidth',1,'MarkerFaceColor',...
        'black','DisplayName','Bar Type 9');
end
if isempty(t6)~=1
    figure(4)
    plot(dispVar6x,dispVar6y,'m o','linewidth',1,'MarkerFaceColor',...
        'magenta','DisplayName','Bar Type 10');
end
if isempty(t7)~=1
    figure(4)
    plot(dispVar7x,dispVar7y,'o','linewidth',1,'MarkerFaceColor',...
        '[0.255 0.069 0]','DisplayName','Bar Type 12');
end
        
%------------------------- end column plot ---------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Not reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagrama1(:,2);
y=diagrama1(:,1);

x2=diagrama2(:,2);
y2=diagrama2(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagrama1(:,4);
yfr=diagrama1(:,3);
npuntos=length(y);

xfr2=diagrama2(:,4);
yfr2=diagrama2(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%% Direction 1-2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=load_conditions(:,3);
ycondicion=load_conditions(:,2);
figure(7)
plot(x,y,'k')
hold on
plot(x2,y2,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on
plot(xfr2,yfr2,'r','DisplayName',' ')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.1,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([-diagrama1(fix(npuntos/2)+1,2)*1.5 diagrama1(fix(npuntos/2)+1,2)*1.5 ...
    diagrama1(1,1) diagrama1(npuntos,1)]);
xlabel('Flexure moment')
ylabel('Axial Force')
title('Interaction diagrams (Rebar) - X-axis')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Not reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=diagrama1(:,6);
y=diagrama1(:,1);

x2=diagrama2(:,6);
y2=diagrama2(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfr=diagrama1(:,8);
yfr=diagrama1(:,7);

xfr2=diagrama2(:,8);
yfr2=diagrama2(:,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcondicion=load_conditions(:,4);
ycondicion=load_conditions(:,2);
figure(8)
plot(x,y,'k')
hold on
plot(x2,y2,'k')
legend('Nominal')
hold on
plot(xfr,yfr,'r','DisplayName','Reduced')
hold on
plot(xfr2,yfr2,'r','DisplayName',' ')
hold on
plot(xcondicion,ycondicion,'r o','linewidth',0.05,'MarkerFaceColor','red',...
    'DisplayName','Load Condition')
axis([-diagrama1(fix(npuntos/2)+1,6)*1.5 diagrama1(fix(npuntos/2)+1,6)*1.5 ...
    diagrama1(1,1) diagrama1(npuntos,1)]);
xlabel('Flexure Bending')
ylabel('Axial Force')
title('Interaction diagrams (Rebar) - Y-axis')
grid on
