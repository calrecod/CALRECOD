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

b=50;%cm
h=50;
E=2.1e6;%kg/cm2
npuntos=50;

numberRebars_hdimension=4;
numberRebars_bdimension=3;


nv=2*numberRebars_hdimension+2*numberRebars_bdimension;
ov=12;
dv=ov/8*2.54; %cm
av=pi/4*dv^2;
RebarTypeIndex=7;

concreteCover=[4 4]; %cm

[dispositionRebar]=RebarDisposition(b,...
    h,concreteCover,dv,nv,numberRebars_hdimension,numberRebars_bdimension);

load_conditions=[1 -30 51 40]; % [num-condition, Pu, Mux, Muy] Ton,m

fdpc=280*0.85;
betac=0.85;

As=nv*av;

[diagrama,maxef,eficiencia,cxy]=diagramasDisposicion(As,b,h,E,npuntos,...
               fdpc,nv,betac,ov,av,dispositionRebar,load_conditions);
 
bestArrangement=zeros(nv,1)+RebarTypeIndex;

diagramsFinalRebarCols(load_conditions,diagrama,dispositionRebar,...
                h,b,bestArrangement);