% BS_CostComparison_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compare the construction cost variation between symmetrical and
%    asymmetrical design for a column according to the variation in the
%    constructability factor fo assembly CFA.
%
%    Notes: Function optrebarColsSymIntSurf is used for the symmetrical
%           rebar design and AssemRebarSymAsymIntSurf is used for the
%           asymmetrical design.
%
%----------------------------------------------------------------
% CREATED:       L.F.Veduzco    2022-06-26
%                Faculty of Engineering
%                Autonomous University of Queretaro
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clc
clear all

%% GEOMETRY
height=400; %cm
fy=4200; % yield stress of rebars 

rec=[4 4]; % [coverx covery] (cm)
E=2.1e6;

%% ADITIONAL PARAMETERS

pu_cols_isr=[28.93]./18;

condition_cracking="Cracked"; % "Cracked" or "Non-cracked"

RebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];

%% Ductility demand
ductility=3;

npdiag=40;

%% Main loop

hglobal=[];

b=40; % cm
h=50; % cm

load_conditions=[1,-11149.6973313160,904768.359104415,2117662.19094098;
                 2,-9601.53706971447,-807897.945266353,-2765903.02800142];
fc=280; % Kg/cm2

fdpc=0.85*fc;
betac=1.05-1400/fc;
if betac>0.85
    betac=0.85;
elseif betac<0.65
    betac=0.65;
end
% --------------------------------------------------------------------
% ISR optimization
% --------------------------------------------------------------------

[b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_01,t_value_03,...
    cxp,bestLoad]=isrColsInterSurf(pu_cols_isr,height,b,h,rec,fy,fc,...
    load_conditions,0.0078,ductility,0,0);

% --------------------------------------------------------------------
% RP-1 (Sym-basic)
% --------------------------------------------------------------------

% ------------------------------------------------------------------
% Symmetrical primary design (Sym-basic)
% ------------------------------------------------------------------
puColBuild=[1.2,0.9,128,214.75,0.04,0.105,7];

[Mr_col,h1,OriginalArea,lowestCost,ovMostEc,nvEc,maxEfEc,...
bestArrangement,bestDisposition,nvxy,bestdiagram,bestLoad]=...
optrebarColsSymIntSurf(b,h,rec,Ac_sec_elem,E,npdiag,fdpc,betac,...
RebarAvailable,0.0078,height,load_conditions,puColBuild,1);

OriginalArea

% --------------------------------------------------------------------
% RP-Assemm (1 pack)
% --------------------------------------------------------------------
dataCFA=[0,1,1,2]; % Data for the computation of the constructability of the
                   % rebar designs. Max and Min values of constructability  
                    % for the rebar designs and weight values of uniformity
                    %  of number of rebars and number of rebar diameters, 
                    % respectively
[h2,bestAreaAssem,bestEf,bestdiagram,bestArrangement,bestDisposition,...
bestMr,bestcxy,bestCP,bestCost,bestLoad,bestCFA2]=AssemRebarSymAsymIntSurf(b,...
h,rec,Ac_sec_elem,E,npdiag,fdpc,betac,height,0.0078,...
RebarAvailable,load_conditions,ductility,puColBuild,dataCFA,0);

bestAreaAssem

% Buildability

% COST AND VOLUME OF REBAR------------------------------------------
% ------------------------------------------------------------------
y1=[OriginalArea,bestAreaAssem].*height*0.0078;

%% Bar diagrams to compare results
x=categorical({'Symmetrical','Asymmetrical'});
figure(10)
bar(x,y1,'black')
hold on
title('Weight of the reinforcing steel in the columns')
xlabel('Type of reinforcement')
ylabel('Weight of the structure (Kgf)')

y2=[lowestCost,bestCost];

figure(11)
bar(x,y2,'black')
hold on
title('Construction cost of the structure: Steel rebar')
xlabel('Type of reinforcement')
ylabel('Cost')

y3=[1,bestCFA2];

figure(12)
bar(x,y3,'black')
hold on
title('Constructability of reinforcing steel rebar design (CFA)')
xlabel('Type of reinforcement')
ylabel('Complexity Factor of Assembly (CFA)')

% ---------------------------------------------------------------------
