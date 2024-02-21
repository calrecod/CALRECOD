% NeutralAxisRecTbeams_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compare the variation of bending resistance between a T-beam
%    cross-section and a rectangular cross-section as the neutral axis
%    depth changes.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2023-03-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% Geometry of the T-beam cross-section
bp=20; % web width (cm) 
ht=40; % total height (cm)
ba=60; % flange width (cm) 
ha=12; % flange height or thickness (cm)
span=500; % cm

cover=4; % lateral concrete cover

%% Materials
fc=250; % Kg/cm2
fy=4200; % Yield stress of steel reinforcement (Kg/cm2)
Es=fy/0.0021; % Modulus of elasticity of the reinforcing steel
fdpc=fc*0.85;

%% Additional general parameters
beta1=0.85; % For the reduction of the neutral axis depth as specified
            % by code

%% Rebar
dispositionRebar=[-7 -17;
                   0 -17;
                   7 -17;
                   -7 17;
                   7 17];
               
rebarType=[4;4;4;4;4];

% Database of the commercially available rebar
rebarAvailable=[4 4/8.*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];     

As=sum(rebarAvailable(rebarType,2).^2.*pi./4); % Total rebar area

%% Mr in function of the neutral axis (T-beam cross-section)
ce=[];
valuesMr=[];
for c=0:0.1:ht
    a=beta1*c;
    eMecVar=eleMecanicosRebarTBeams(c,a,fdpc,ha,ba,bp,ht,span,Es,rebarType,...
                                    dispositionRebar,rebarAvailable);
                                
    Mr=eMecVar(1,2)+eMecVar(2,2);
    
    valuesMr=[valuesMr,Mr];
    ce=[ce,c];
end

figure(2)
plot(ce,valuesMr,'k -')
legend('T-beam section')
xlabel('Neutral axis depth c (cm)')
ylabel('Bending resistance (Kg-cm)')
title('Bending resistance in function of the neutral axis depth for a beam cross-section')
grid on
hold on

%% Geometry of a rectangular beam cross-section
b=20; % section's width (cm) 
h=40; % section's height (cm)

%% Rebar
rebarTypeT=[4;4;4];
rebarTypeC=[4;4];
  
%% Mr in function of the neutral axis (rectangular cross-section)
ce=[];
valuesMr=[];
for c=0:0.1:h
    a=beta1*c;
    
    eMecVar=eleMecanicosRebarBeams(c,a,fdpc,h,b,cover,Es,rebarTypeT,...
                            rebarTypeC,dispositionRebar,rebarAvailable);
                        
    Mr=eMecVar(1,2)+eMecVar(2,2); 
    
    valuesMr=[valuesMr,Mr];
    ce=[ce,c];
end

figure(2)
plot(ce,valuesMr,'b -','DisplayName','Rectangular section')
grid on
hold on