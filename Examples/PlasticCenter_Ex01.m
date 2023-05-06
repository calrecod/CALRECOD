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

%% Column's geometry
b=40; % cross-section width dimension (cm)
h=40; % cross-section height dimension (cm)

concreteCover=[4 4]; % concrete cover in both cross-section directions (cm)

%% Materials
fc=300; % concrete's compressive strength (Kg/cm2)
fdpc=fc*0.85; % reduced concrete's strength (Kg/cm2)
fy=4200; % yield stress of the reinforcing steel (Kg/cm2)
           
%% Rebar data
% Commercially available rebar diameters. [eight-of-an-inch, diam(cm)]
rebarAvailable=[4 4/8*2.54;
                    5 5/8*2.54;
                    6 6/8*2.54;
                    8 8/8*2.54;
                    9 9/8*2.54;
                    10 10/8*2.54;
                    12 12/8*2.54];

% number of rebars at each cross-section's boundary
numberRebars1=9;
numberRebars2=3;
numberRebars3=5;
numberRebars4=2;

% rebar diameter at each cross-section's boundary
RebarTypeIndex1=4;
RebarTypeIndex2=3;
RebarTypeIndex3=5;
RebarTypeIndex4=3;
                
% Total number of rebars placed over the cross-section
nv=numberRebars1+numberRebars2+numberRebars3+numberRebars4; % number of 
                                                            % rebars 
% Rebar distribution
[dispositionRebar,separacion_hor1,separacion_hor2,...
separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
h,concreteCover,nv,numberRebars1,numberRebars2,...
numberRebars3,numberRebars4,rebarAvailable,RebarTypeIndex1,...
RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4);

%% Computation of the Plastic Center with respect to the X-axis
rebarTypeslist(1:numberRebars1)=RebarTypeIndex1;
rebarTypeslist(numberRebars1+1:numberRebars1+numberRebars2)=RebarTypeIndex2;
rebarTypeslist(numberRebars1+numberRebars2+1:numberRebars1+numberRebars2+...
    numberRebars3)=RebarTypeIndex3;
rebarTypeslist(numberRebars1+numberRebars2+numberRebars3+1:nv)=RebarTypeIndex4;

[PCX]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)
                        
%% Computation of the Plastic Center with respect to the Y-axis
rebar=dispositionRebar; 

% Invert rebar local coordinates (the cross-section is rotated 90°)
dispositionRebar(:,1)=-rebar(:,2);
dispositionRebar(:,2)=rebar(:,1);

% Invert cross-section dimensions
dimensionesColumna=[b h];
h=dimensionesColumna(1);
b=dimensionesColumna(2);  

% Invert rebar diameters
rebarTypeslist(1:numberRebars3)=RebarTypeIndex3;
rebarTypeslist(numberRebars3+1:numberRebars3+numberRebars4)=RebarTypeIndex4;
rebarTypeslist(numberRebars3+numberRebars4+1:numberRebars3+numberRebars4+...
    numberRebars2)=RebarTypeIndex2;
rebarTypeslist(numberRebars3+numberRebars4+numberRebars2+1:nv)=RebarTypeIndex1;

[PCY]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,rebarTypeslist,...
                            rebarAvailable)   
                        
%% Plotting results
% Rebar diameters' indices
barTypes1(1:numberRebars1)=RebarTypeIndex1;
barTypes1(numberRebars1+1:numberRebars1+numberRebars2)=RebarTypeIndex2;
barTypes2(1:numberRebars3)=RebarTypeIndex3;
barTypes2(numberRebars3+1:numberRebars3+numberRebars4)=RebarTypeIndex4;

% Calling function "beamReinforcedSection"
beamReinforcedSection(h,b,rebar,barTypes1,barTypes2)

% Coordinates of the Plastic Center
xpc=h/2-PCX;
ypc=b/2-PCY;

% Ploting the Plastic Center location
PCXText=num2str(xpc); % To convert the numerical values in strings
PCYText=num2str(ypc);
PCXY=strcat('( ',PCXText,', ',PCYText,' )'); % To concatenate strings
figure(2)
plot(ypc+2.5,xpc-2.5,'o','MarkerFaceColor','magenta','DisplayName','PC')
text(ypc,xpc,PCXY) % To plot the coordinates as text
