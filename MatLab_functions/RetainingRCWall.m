function [compliedRestric,areaWall,linearWeigthWall,tippingFS,slideFS,...
    LCap_FS,sepheel,efHeel,sepfoot,effoot,septrunk,eftrunk]=RetainingRCWall...
    (foot,heel,hf,b,FiFill,H,D,m1,m2,wvFill,beta,FiBackFill,alfa,FiFound,...
    fc,fy,wvc,ductility,qadm,minFSqadm,SlideSF,TippingSF,typeRebar,sepMinRebars,...
    maxEf,qaf,qab,qs,LF_DL)
        
%-------------------------------------------------------------------------
% Syntax:
% [compliedRestric,areaWall,linearWeigthWall,tippingFS,slideFS,...
%  LCap_FS,sepheel,efHeel,sepfoot,effoot,septrunk,eftrunk]=RetainingRCWall...
%  (foot,heel,hf,b,FiFill,H,D,m1,m2,wvFill,beta,FiBackFill,alfa,FiFound,...
%  fc,fy,wvc,ductility,qadm,minFSqadm,SlideSF,TippingSF,typeRebar,sepMinRebars,...
%  maxEf,qaf,qab,qs,LF_DL)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%-------------------------------------------------------------------------
% PURPOSE: To analyse and design the reinforcement of a concrete retaining 
% wall for certain given cross-section dimensions, design Safety Factors
% and an admisible soil's bearing load capacity as well as other soil
% fill's mechanical properties.
% 
% OUTPUT: areaWall:             is the wall's cross-sectional area
%
%         linearWeightWall:     is the linear weight (per cm) of the
%                               designed retaining wall
%
%         tippingFS:            is the final tipping Safety Factor for the
%                               designed wall
%
%         slideFS:              is the final slide Safety Factor for the
%                               designed wall
%
%         LCap_FS:              is the final Safety Factor against the
%                               soil's bearing load capacity for the
%                               designed wall
%
%         sepheel:              is the final rebar separation for the
%                               designed reinforcement in the concrete
%                               wall's heel
%
%         efHeel:               is the final structural efficiency for the
%                               designed wall's heel
%
%         sepfoot:              is the final rebar separation for the
%                               designed reinforcement in the concrete
%                               wall's foot
%
%         effoot:               is the final structural efficiency for the
%                               designed wall's foot
%
%         septrunk:             is the final rebar separation for the
%                               designed reinforcement in the concrete
%                               wall's trunk
%
%         eftrunk:              is the final structural efficiency for the
%                               designed wall's trunk
%
% INPUT:  foot:                 is the cross-sectional length of the wall's
%                               foot
%
%         heel:                 is the cross-sectional length dimension of
%                               the wall's heel
%
%         hf:                   is the cross-sectional width of the wall's 
%                               heel and foot
%
%         b:                    is the upper cross-section wall's trunk
%                               width
%
%         FiFill:               is the fill soil's friction angle
%
%         H:                    is the wall's stem height dimension
%
%         D:                    is the soil's bakc fill depth
%
%         m1, m2:               are the wall's dowels' from and back slopes
%
%         wvFill:               is the soil fill's volumetric weight
%
%         beta:                 is the front fill soil's upper-grade angle
%
%         FiBackFill:           is the back fill soil's friction angle
%
%         alfa:                 is the back fill soil's upper-grade angle
%
%         FiFound:              is the foundation soil's friction angle
%
%         fc:                   is the concrete compressive strength f'c
%                               (in units Kg/cm^2)
%
%         wvc:                  is the concrete volumetric weight (Kg/cm^3)
%
%         ductility:            is the level of ductility demand for the
%                               design of the reinforcement
%
%         qadm:                 is the soil's bearing load capacity, in 
%                               units (Kg/cm^2)
%
%         minFSqadm:            is the design Safety Factor against the
%                               soil's bearing load capacity
%
%         SlideSF:              is the design Safety Factor against the
%                               slide forces over the wall
%
%         TippingSF:            is the design Safety Factor against the
%                               tipping forces over the wall
% 
%         typeRebar:            is the rebar eigth-of-an-inch to be used
%                               for the wall's reinforcement
%
%         sepMinRebars:         is the minimum rebar separation restriction
%                               for the wall's reinforcement
%
%         maxEf:                is the critical structural efficiency for 
%                               the design of the wall's heel, foot and 
%                               trunk
%
%         qaf, qab:             is the front surcharge and the back 
%                               surcharge, respectively (in units Kg/cm^2) 
%
%         qs:                   is the linear load that the wall may
%                               support over the top along its length
%
%         LF_DL:                is the Dead Load Design Factor
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%% Additional important dimensions
tb=b+1/m1*H+1/m2*H;
B=foot+tb+heel;

%% Soil fill
KaFill=tan(deg2rad(45-FiFill*0.5))^2;
KpFill=tan(deg2rad(45+FiFill*0.5))^2;

deltaFill=2/3*FiFill;

Ea=wvFill/2*H^2*KaFill;
KoFill=1-sin(deg2rad(FiFill));

%% Backing ground fill
KaBackFill=tan(deg2rad(45-FiBackFill*0.5))^2;
KpBackFill=tan(deg2rad(45+FiBackFill*0.5))^2;

deltaBackFill=2/3*FiBackFill;
wvBackFill=wvFill;
Ep=wvBackFill/2*D^2*KaBackFill;
KoBackFill=1-sin(deg2rad(FiBackFill));

%% Foundation's soil
deltaFound=2/3*FiFound;

%% Concrete and reinforcement
if fc<=280 
    beta1=0.85;
elseif fc>280
    beta1=1.05-fc/1400;
    if beta1<0.65
        beta1=0.65;
    end
end  
fdpc=0.85*fc;

% Rebar separation
ab=(typeRebar/8*2.54)^2*pi/4;
as1=660*(hf-5+5)/(fy*(hf+100));
if ab/as1>(hf*3.5) || ab/as1>50
    if hf*3.5>50
        
        sepMaxRebars=50;
    else
        sepMaxRebars=fix(hf*3.5);
    end
else
    sepMaxRebars=fix(ab/as1);
end

%% Table of analysis
Me=[]; Md=[];

% Concrete wedges ------------------------------------------
% wedge1                            % wedge2 
area1=0.5*1/m1*H*H;                 area2=b*H;
linearWeight1=area1*wvc;            linearWeight2=area2*wvc;
x1=foot+2/3*(1/m1*H);               x2=foot+(1/m1*H)+0.5*b;

% Stabilizing moment 1              % Stabilizing moment 2
Me1=linearWeight1*x1;               Me2=linearWeight2*x2;
Me=[Me;Me1];                        Me=[Me;Me2];

% ----------------------------------------------------------
% wedge3                            % wedge4
area3=0.5*1/m2*H*H;                 area4=hf*B;
linearWeight3=area3*wvc;            linearWeight4=area4*wvc;
x3=foot+(1/m1*H)+b+1/3*(1/m2*H);    x4=0.5*B;

% Stabilizing moment 3              % Stabilizing moment 4
Me3=linearWeight3*x3;               Me4=linearWeight4*x4;
Me=[Me;Me3];                        Me=[Me;Me4];

% Soil fill wedges -----------------------------------------
% wedge5                            % wedge 6
area5=foot*D;                       area6=0.5*1/m1*D*D;
linearWeight5=area5*wvBackFill;     linearWeight6=area6*wvBackFill;
x5=0.5*foot;                        x6=foot+1/3*(1/m1*D);

% Stabilizing moment 5              % Stabilizing moment 6
Me5=linearWeight5*x5;               Me6=linearWeight6*x6;
Me=[Me;Me5];                        Me=[Me;Me6];

% ----------------------------------------------------------
% wedge 7                           % wedge 8
area7=0.5*1/m2*H*H;                 area8=heel*H;
linearWeight7=area7*wvFill;         linearWeight8=area8*wvFill;
x7=foot+(1/m1*H)+b+2/3*(1/m2*H);    x8=foot+tb+0.5*heel;

% Stabilizing moment 7              % Stabilizing moment 8
Me7=linearWeight7*x7;               Me8=linearWeight8*x8;
Me=[Me;Me7];                        Me=[Me;Me8];

% ----------------------------------------------------------
% Eay                               % Epx
Eay=Ea*sin(deg2rad(beta));          Epx=Ep*cos(deg2rad(alfa));
xeay=foot+tb+2/3*(1/m2*H);          yepx=hf+1/3*D;

% Stabilizing moment 9              % Stabilizing moment 10
Meay=Eay*xeay;                      Mepx=Epx*yepx;
Me=[Me;Meay];                       Me=[Me;Mepx];

% ----------------------------------------------------------
% Epy                               % qa (back)
Epy=Ep*sin(deg2rad(alfa));          qb=qab*D*KoFill;
xepy=foot+1/3*(1/m1*D);             yqb=hf+0.5*D;

% Stabilizing moment 11             % Stabilizing moment 12
Mepy=Epy*xepy;                      Meqb=qb*yqb;
Me=[Me;Mepy];                       Me=[Me;Meqb];

% ----------------------------------------------------------
% Eax                               % qa (front)
Eax=Ea*cos(deg2rad(beta));          qf=qaf*H*KoFill;
yeax=hf+1/3*H;                      yqf=hf+0.5*H;

% Destabilizing moment 1            % Destabilizing moment 2
Md1=Eax*yeax;                       Md2=qf*yqf;
Md=[Md;Md1];                        Md=[Md;Md2];

%% Transversal area of wall - Objective function
areaWall=area1+area2+area3+area4;
linearWeigthWall=areaWall*wvc;

%% Results of analysis
% Tipping Safety Factor -----------------------------------------------
tippingFS=0.7*sum(Me)/(LF_DL*sum(Md));
if tippingFS<TippingSF
    tippingRestriction=0;
    % disp('The design does not comply by Tipping Safety Factor')
else
    tippingRestriction=1;
end
% Slide Safety Factor -------------------------------------------------
Nqaf=qaf*(heel+(1/m2)*H);

% Normal vertical force that contributes to 
% the sliding resistance by soil friction 
% (Live Loads are not considered)
N=linearWeight1+linearWeight2+linearWeight3+linearWeight4+linearWeight5+...
    linearWeight6+linearWeight7+linearWeight8+Eay+Nqaf;

Fr=N*tan(deg2rad(deltaFound));

Fd=Eax+qf;

slideFS=0.9*(Fr)/(LF_DL*Fd);

if slideFS<SlideSF
    % disp('The design does not comply by Sliding Safety Factor')
    slideRestriction=0;
else
    slideRestriction=1;
end

% Bearing load capacity - Contact pressures ---------------------------
% Service Loads are considered
xR=(sum(Me)-sum(Md))/N;
ecc=B/2-xR;

% Normal vertical force that the soild widthstands 
% (Live Loads and Dead Loads are considered - Service Load Case)
Nqab=qab*(foot+(1/m1)*D);

Nq=linearWeight1+linearWeight2+linearWeight3+linearWeight4+linearWeight5+...
    linearWeight6+linearWeight7+linearWeight8+qs+Eay+Epy+Nqaf+Nqab;

sigmaLeft=(Nq)/B*(1+(6*ecc/B));
sigmaRight=(Nq)/B*(1-(6*ecc/B));

LCap_FS=qadm/max([sigmaLeft,sigmaRight]);
if LCap_FS<minFSqadm
    % disp('The soil load capacity has been overpast')
    qadmRestriction=0;
else
    qadmRestriction=1;
end

%% Mechanic Elements
% Heel
Qheel=LF_DL*min([sigmaRight,sigmaLeft])*heel+...
       LF_DL*abs(sigmaLeft-sigmaRight)/B*heel*0.5*heel;
   
Vheel=Qheel/hf;
Mheel=LF_DL*min([sigmaRight,sigmaLeft])*0.5*heel^2+...
       LF_DL*abs(sigmaLeft-sigmaRight)/B/6*heel^3;

% Foot
Qfoot=LF_DL*min([sigmaRight,sigmaLeft])*(heel+tb)+...
       LF_DL*abs(sigmaLeft-sigmaRight)/B*(heel+tb)^2*0.5;
   
Vfoot=Qfoot/hf;
Mfoot=LF_DL*min([sigmaRight,sigmaLeft])*0.5*(heel+tb)^2+...
       LF_DL*abs(sigmaLeft-sigmaRight)/B/6*(heel+tb)^3;
   
% Trunk
Qtrunk=LF_DL*(Eax+qf);
   
Vtrunk=Qtrunk/tb;
Mtrunk=LF_DL*(Eax*H/3+qf*H/2);

%% Design of reinforcement
pmin=0.0026;
if ductility==1 || ductility==2
    pmax=(0.9*fdpc/fy*(6000*beta1)/(fy+6000));
elseif ductility==3
    pmax=(0.75*fdpc/fy*(6000*beta1)/(fy+6000));
end

% Heel ---------------------------------------------------------------

qheel=(1-sqrt(1-4*0.5*Mheel/((hf-5)^2*0.9*fdpc)))/(2*0.5);
pheel=qheel*fdpc/fy;

if pheel<pmin
    pheel=pmin;
    amaxConditionHeel=1;
elseif pheel>pmax
    pheel=pmax;
    % disp('The heel dimension is too small')
    amaxConditionHeel=0; % rebar area is not enough - dimensions too small

elseif pheel<pmax && pheel>pmin
    amaxConditionHeel=1;
end

qheel=pheel*fy/fdpc;

Asheel=pheel*100*(hf-5);
nbheel=ceil(Asheel/ab);
if nbheel<2
    nbheel=2;
end
sepheel=(100-nbheel*(typeRebar/8*2.54))/(nbheel-1)-...
         mod((100-nbheel*(typeRebar/8*2.54))/(nbheel-1),5);
     
if sepheel<sepMinRebars % separation is not enough
    % disp('The rebar separation in the heel is too small')
    sepConditionHeel=0;
elseif sepheel>sepMaxRebars 
    sepheel=sepMaxRebars;
    sepConditionHeel=1;
    
elseif sepheel<=sepMaxRebars && sepheel>=sepMinRebars
    sepConditionHeel=1;
end

Mrheel=0.9*(pheel*(hf-5)^2*fy*(1-0.5*qheel));
efHeel=Mheel/Mrheel;

if efHeel>maxEf
    % disp('The resistance of the heel by flexure is not enough')
    efConditionHeel=0; % resistance is not enough
else
    efConditionHeel=1;
end

% Foot ---------------------------------------------------------------

qfoot=(1-sqrt(1-4*0.5*Mfoot/((hf-5)^2*0.9*fdpc)))/(2*0.5);
pfoot=qfoot*fdpc/fy;

if pfoot<pmin
    pfoot=pmin;
    amaxConditionFoot=1;
elseif pfoot>pmax
    pfoot=pmax;
    % disp('The foot dimension is too small')
    amaxConditionFoot=0; % rebar area is not enough - dimensions too small

elseif pfoot<pmax && pfoot>pmin
    amaxConditionFoot=1;
end

qfoot=pfoot*fy/fdpc;

Asfoot=pfoot*100*(hf-5);

nbfoot=ceil(Asfoot/ab);
if nbfoot<2
    nbfoot=2;
end  
sepfoot=(100-nbfoot*(typeRebar/8*2.54))/(nbfoot-1)-...
         mod((100-nbfoot*(typeRebar/8*2.54))/(nbfoot-1),5);
     
if sepfoot<sepMinRebars % separation is not enough
    % disp('The rebar separation in the foot is too small')
    sepConditionFoot=0;
elseif sepfoot>sepMaxRebars 
    sepfoot=sepMaxRebars;
    sepConditionFoot=1;
elseif sepfoot<=sepMaxRebars && sepfoot>=sepMinRebars
    sepConditionFoot=1;
end

Mrfoot=0.9*(pfoot*(hf-5)^2*fy*(1-0.5*qfoot));
effoot=Mfoot/Mrfoot;

if effoot>maxEf
    % disp('The resistance of the foot under flexure is not enough')
    efConditionFoot=0; % resistance is not enough
else
    efConditionFoot=1;
end

% Trunk --------------------------------------------------------------

qtrunk=(1-sqrt(1-4*0.5*Mtrunk/((tb-5)^2*0.9*fdpc)))/(2*0.5);
ptrunk=qtrunk*fdpc/fy;

if ptrunk<pmin
    ptrunk=pmin;
    amaxConditionTrunk=1;
elseif ptrunk>pmax
    ptrunk=pmax;
    % disp('The trunk width dimension at the base is too small')
    amaxConditionTrunk=0; % rebar area is not enough - dimensions too small

elseif ptrunk<pmax && ptrunk>pmin
    amaxConditionTrunk=1;
end

qtrunk=ptrunk*fy/fdpc;

Astrunk=ptrunk*100*(tb-5);
nbtrunk=ceil(Astrunk/ab);
if nbtrunk<2
    nbtrunk=2;
end  

septrunk=(100-nbtrunk*(typeRebar/8*2.54))/(nbtrunk-1)-...
         mod((100-nbtrunk*(typeRebar/8*2.54))/(nbtrunk-1),5);
   
if septrunk<sepMinRebars % separation is not enough
    sepConditionTrunk=0;
    % disp('The rebar separation in the trunk is too small')
elseif septrunk>sepMaxRebars 
    septrunk=sepMaxRebars;
    sepConditionTrunk=1;
elseif septrunk<=sepMaxRebars && septrunk>=sepMinRebars
    sepConditionTrunk=1;
end

Mrtrunk=0.9*(ptrunk*(tb-5)^2*fy*(1-0.5*qtrunk));
eftrunk=Mtrunk/Mrtrunk;

if eftrunk>maxEf
    % disp('The resistance of the trunk under flexure is not enough')
    efConditionTrunk=0; % resistance is not enough
else
    efConditionTrunk=1;
end

restrictions=[tippingRestriction,slideRestriction,qadmRestriction,...
        amaxConditionHeel,sepConditionHeel,efConditionHeel,amaxConditionFoot,...
        sepConditionFoot,efConditionFoot,amaxConditionTrunk,sepConditionTrunk,...
        efConditionTrunk];

if sum(restrictions)==12
    % disp('The design complies with all design restrictions')
    compliedRestric=1;
else
    compliedRestric=0;
end
