% StabilitySlenderRCols_Rebar_Ex02
%----------------------------------------------------------------
% PURPOSE 
%    To analyse the variation in structural efficiency of an symmetrically
%    reinforced rectangular cross-section when considering 
%    cracking-cross-section mechanisms and Geometrical Non-linearities
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-04-16
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Geometry
b=40;
h=40;
dimensionsColumn=[b h];
height=500; % column's length

%% Materials
fc=300; % concrete's compressive strength
Ec=14000*sqrt(fc); % Modulus of elasticity of the concrete

fdpc=fc*0.85;
betac=0.85;

fy=4200; % yield stress of the reinforcing steel (Kg/cm2)
E=2.1e6; % Modulus of elasticity of the reinforcing steel (kg/cm2)

%% Additional structural parameters
k=2; % slenderness factor
npdiag=50; % number of points to be computed for the int. diagrams
concreteCover=[4 4]; % cm

%% Loads
load_conditions=[1 -15000 28e5 22e5]; % [nload, Pu, Mx, My]

P=load_conditions(2);
Mx=load_conditions(3);
My=load_conditions(4);
Vx=0;
Vy=0;

%% Rebar data
% Commmercially available rebar diameters (eight-of-an-inch)
rebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];

% Distribution of rebar
numberRebars_hdimension=4;
numberRebars_bdimension=5;
nv=2*numberRebars_hdimension+2*numberRebars_bdimension;

% Rebar diameter and area
RebarTypeIndex=4;

ov=rebarAvailable(RebarTypeIndex,1);
dv=rebarAvailable(RebarTypeIndex,2); % cm
av=pi/4*dv^2;

% Total rebar area over the cross-section
As=nv*av;

ab1=numberRebars_bdimension*rebarAvailable(RebarTypeIndex,2)^2*pi/4;
ab2=numberRebars_hdimension*rebarAvailable(RebarTypeIndex,2)^2*pi/4;

% Rebar coordinates over the cross-section
[dispositionRebar]=RebarDisposition(b,h,concreteCover,dv,nv,...
    numberRebars_hdimension,numberRebars_bdimension);

%% Compute the column's interaction diagram
[diagrama,cPoints,Poc,Pot]=diagramRColumnSymRebar(As,b,h,E,npdiag,...
                       fdpc,nv,betac,ov,av,dispositionRebar);
                   
%% Column's structural resistance efficiency 
% according to the initially applied load combination
[maxef01,eficiencia01,cxy01]=effRecColsLinearSearch...
        (diagrama,load_conditions,Pot,Poc,cPoints)     
    
bestArrangement=zeros(nv,1)+RebarTypeIndex;

%% Plot the interaction diagram for better assessment 
diagramsFinalRebarCols(load_conditions,diagrama,dispositionRebar,...
                        h,b,bestArrangement);
                    
%% Modified inertia momentum of the cross-section
% considering mechanism with the transformed cross-section method
eccentricityXY=[abs(Mx/P),abs(My/P)]; % cm
conditionCrack="Cracked";

t1bar=ab1/(b-2*concreteCover(1));
t2bar=ab2/(b-2*concreteCover(1));
InertiaXY=[1/12*b*h^2,1/12*h*b^3];
[InertiaXY,Atransf_xy,elimxy]=CrackingColumnsSym(h,b,fdpc,concreteCover,...
          t1bar,eccentricityXY,t2bar,P,cxy01,conditionCrack,E);
 
% Determine the Euler's critical load: the smaller inertia value of both
% axis directions must be taken as reference
if InertiaXY(1)<=InertiaXY(2)
    I=InertiaXY(1);
else
    I=InertiaXY(2);
end

Pcr=pi^2*Ec*I/(k*height)^2;

%% Amplification bending moment for each axis direction
[Deltax,Mampx]=MomAmpColsGeomNL(fc,k,InertiaXY(1),height,Vx,P,Mx,b,h,0); % In the x-axis
h=dimensionsColumn(1);
b=dimensionsColumn(2);   
[Deltay,Mampy]=MomAmpColsGeomNL(fc,k,InertiaXY(2),height,Vy,P,My,b,h,0); % In the y-axis

%% Amplified load bending moments
load_conditions_modif=[1 P Mampx Mampy];

%% Reduction in resistance due to the amplified moment by
% geometrical non-linearities (P-Delta effects)
[maxef02,eficiencia02,cxy02]=effRecColsLinearSearch...
        (diagrama,load_conditions_modif,Pot,Poc,cPoints)

%% Plot the interaction diagram for better assessment of amplified moments 
diagramsFinalRebarCols(load_conditions_modif,diagrama,dispositionRebar,...
                        h,b,bestArrangement);