% EfficiencyRebarColumnDesign_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the structural resistance efficiency of a symmetrical
%    rebar design of a rectangular column subject to biaxial bending
%
%    Note: function diagramasDisposicion is the only one required to
%          determine such structural efficiency and function 
%          diagramsFinalRebarCols is used only to plot the interaction
%          diagrams and the rebar cross-section design
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clear all
clc

%% Geometry
b=30; % cross-section width (cm)
h=25; % cross-section height (cm)

%% Materials
E=2.1e6; % Modulus of Elasticity of the reinforcing steel (kg/cm2)
fc=250; % Concrete's compressive strength
fdpc=fc*0.85; % reduced f'c
betac=0.85;

%% Additional parameters
npdiag=30; % Number of points to be computed for the interaction diagram
concreteCover=[3 3]; %cm

%% Loads
load_conditions=[1 -9.95e3 5.26e5 0.97e5]; % [nload, Pu, Mx, My]

%% Rebar data
% Commerically available rebar diameters
rebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];
                
% Number of rebars placed vertically and horizontally
numberRebars_hdimension=0; % number of rebars placed vertically (per side)
numberRebars_bdimension=3; % number of rebars placed horizontally

% Total number of rebars placed over the cross-section
nv=2*numberRebars_hdimension+2*numberRebars_bdimension;

RebarTypeIndex=2; % rebar diameter index

ov=rebarAvailable(RebarTypeIndex,1); % eight-of-an-inch rebar
dv=rebarAvailable(RebarTypeIndex,2); % rebar diameter (cm)
av=pi/4*dv^2; % rebar area

% Rebar distribution over the cross-section
[dispositionRebar]=RebarDisposition(b,h,concreteCover,dv,nv,...
    numberRebars_hdimension,numberRebars_bdimension);

%% Additional output: design information
As=nv*av;
disp('Rebar area proposed:')
disp(As)

Asmin=0.01*b*h;
disp('Min rebar area:')
disp(Asmin)

Asmax=0.04*b*h;
disp('Max rebar area:')
disp(Asmax)

%% Structural efficiency
[diagrama,maxef,eficiencia,cxy]=diagramasDisposicion(As,b,h,E,npdiag,...
               fdpc,nv,betac,ov,av,dispositionRebar,load_conditions);

%% Plotting results
bestArrangement=zeros(nv,1)+RebarTypeIndex;

diagramsFinalRebarCols(load_conditions,diagrama,dispositionRebar,...
                h,b,bestArrangement);