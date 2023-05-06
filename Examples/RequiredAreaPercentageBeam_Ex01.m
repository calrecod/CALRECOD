% RequiredAreaPercentageBeam_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the increment in bending moment resistance of a
%    rectangular beam cross-section as the reinforcement in tension and/or
%    compression increases.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-26
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-04-16
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=30;
h=70;

b_rec=4;
h_rec=4;
d=h-h_rec;

%% Materials
fc=250;
factor_fc=0.85;
E=2e6;
fy=4200;

%% Rebar data
disposition_rebar=[-12 42;
                    -6 42;
                    0 42;
                    6 42;
                    12 42;
                    -12 -42;
                    -6 -42;
                    0 -42;
                    6 -42;
                    12 -42];
                
% Available commercial rebar diameters to choose from
%              type  diam     area     
varDisponibles=[4 4/8.*2.54 (4/8*2.54).^2*pi./4;
                5 5/8*2.54 (5/8*2.54).^2*pi./4;
                6 6/8*2.54 (6/8*2.54).^2*pi./4;
                8 8/8*2.54 (8/8*2.54).^2*pi./4;
                9 9/8*2.54 (9/8*2.54).^2*pi./4;
                10 10/8*2.54 (10/8*2.54).^2*pi./4;
                12 12/8*2.54 (12/8*2.54).^2*pi./4];

ast=5*varDisponibles(2,3); % rebar area in tesnion
asc=ast; % rebar area in compression

%% Additional parameters
duct=30; % ductility demand
mmu=3.5; % number of times to increment the initial bending resistance

amax=(0.025*b*d); % max rebar area allowed in tension

%% Main process

% Initial resistance with the current reinforcement quantity
[eftk1,mrast1,ck1]=Efrec2tBeams([1 20e5],fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,0.85);
mrast=mrast1;
mmax=mmu*mrast;    

% Loop
ast_vector=[];
Mrast_vector=[]; 
mrast=0;
while mrast<mmax
    [eftk1,mrast,ck1]=Efrec2tBeams([1 20e5],fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,0.85);
    ast=ast+0.1;
    atotal=ast;
    ast_vector=[ast_vector,atotal];
    Mrast_vector=[Mrast_vector,mrast]; 
end

%% Plotting results
figure(4)
plot(ast_vector,...
   Mrast_vector,'b -','MarkerFaceColor','blue','LineWidth',2.5)
hold on
xlabel('Area of reinforcing steel in tension (cm)')
ylabel('Resistance (Kg-cm)')
title('Resistance of a beam as the area in tension increases')
hold on

