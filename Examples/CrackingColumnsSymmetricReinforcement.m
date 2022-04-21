
clc
clear all

h=50;
b=40;
eccentricityXY=[120,115];
rec=[4,4];
tx=0.005;
ty=0.002;
Pu=28;
cxy=[20,13];
conditionCrack="Cracked";
E=2.1e6;
fdpc=0.85*300;

Inertia0=[b*h^3/12,h*b^3/12]

[InertiaXYmodif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,...
          tx,eccentricityXY,ty,Pu,cxy,conditionCrack,E)