% PlasticCenter_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To determine the location of the Plastic Center of an asymmetrically 
%    reinforced column's cross-section (with respect to the upper outter 
%    concrete fibre) 
%
%    Note: function dispositionRebarAsymmetric is used to determine the
%          rebar distribution as local coordinates according to the given 
%          rebar data
%
%          function PlastiCenterAxis is the one used to compute the plastic
%          center for each cross-section's axis direction (X,Y) separately
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2022-07-27
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clear all
clc

b=40;
h=40;
fc=300;

rebarAvailable=[4 4/8*2.54 0.994;
                    5 5/8*2.54 1.552;
                    6 6/8*2.54 2.235;
                    8 8/8*2.54 3.973;
                    9 9/8*2.54 5.033;
                    10 10/8*2.54 6.225;
                    12 12/8*2.54 8.938];
                
numberRebars1=9;
numberRebars2=2;
numberRebars3=5;
numberRebars4=5;

RebarTypeIndex1=3;
RebarTypeIndex2=3;
RebarTypeIndex3=3;
RebarTypeIndex4=3;
                
dbmin=max([rebarAvailable(RebarTypeIndex1,2),rebarAvailable(RebarTypeIndex2,2),...
           rebarAvailable(RebarTypeIndex3,2),rebarAvailable(RebarTypeIndex4,2)]);
                  
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4; % number of 
                                                            % rebars
concreteCover=[4 4]; %cm
sepMin=max([1.5*2.54 1.5*dbmin]); %cm
                
[dispositionRebar,separacion_hor1,separacion_hor2,...
            separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
            h,sepMin,concreteCover,nv,numberRebars1,numberRebars2,...
            numberRebars3,numberRebars4,rebarAvailable,RebarTypeIndex1,...
            RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4);
      
fdpc=fc*0.85; % (Kg/cm2)
fy=4200; % Kg/cm2

rebarTypeslist=zeros(nv,1);
rebarTypeslist(1:numberRebars1)=RebarTypeIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarTypeIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarTypeIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarTypeIndex4;

% Computation of the Plastic Center with respect to the X-axis
[PCX]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)
                        
% Computation of the Plastic Center with respect to the Y-axis
rebar=[dispositionRebar(:,1) dispositionRebar(:,2)]; 

% Invert rebar local coordinates (the cross-section is rotated 90°)
dispositionRebar(:,1)=-rebar(:,2);
dispositionRebar(:,2)=rebar(:,1);

% Invert cross-section dimensions
dimensionesColumna=[b h];
h=dimensionesColumna(1);
b=dimensionesColumna(2);                   
[PCY]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)            