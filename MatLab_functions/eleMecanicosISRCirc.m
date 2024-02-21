function eleMec=eleMecanicosISRCirc(c,a,fdpc,diam,rec,Es,t)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicosISRCirc(c,a,fdpc,diam,rec,E,t)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the axial load and bending resistance of a 
% reinforced column of circular cross section with an ISR.
% 
% OUTPUT: eleMec:             array containing the sum of resistance forces
%                             (axial and being) as [ sum Fs,   sum Ms;
%                                                   sum Fconc, sum Mconc]
%
% INPUT:  diam:               cross-section diameter
%
%         t:                  is the ISR width
%
%         rec:                is the concrete cover
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code
%
%         Es:                 is the Modulus of Elasticity of the
%                             reinforcing steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

ndt=40; % number of discrete ISRs
teta=360/ndt;
eMecProfile=zeros(ndt,8); % [dA, xposition, yposition, d, eps, eE, Fr, Mr]

% calculate dA for each t area
dA=(pi*(diam-2*rec))*t/ndt;

% localtion of each discrete dt over the section
for i=1:ndt
    eMecProfile(i,1)=dA;
    x=(diam-2*rec)*0.5*cos(deg2rad(teta*0.5+teta*(i-1)));
    eMecProfile(i,2)=x;

    if i<=0.5*ndt
        y=sqrt((diam*0.5)^2-x^2);
        eMecProfile(i,3)=y;
    elseif i>0.5*ndt
        y=-sqrt((diam*0.5)^2-x^2);
        eMecProfile(i,3)=y;
    end
end

sumaM=0;
sumaF=0;
for i=1:ndt
    eMecProfile(i,4)=0.5*diam-eMecProfile(i,3); % momentum distance
    eMecProfile(i,5)=0.003/c*(eMecProfile(i,4)-c);
    
    if (eMecProfile(i,5)<-0.0021)
        eMecProfile(i,5)=-0.0021;
    elseif(eMecProfile(i,5)>0.0021)
        eMecProfile(i,5)=0.0021;
    end
    eMecProfile(i,6)=eMecProfile(i,5)*Es;
    eMecProfile(i,7)=eMecProfile(i,6)*eMecProfile(i,1);
    eMecProfile(i,8)=eMecProfile(i,7)*(eMecProfile(i,4)-0.5*diam);
    
    sumaF=sumaF+eMecProfile(i,7);
    sumaM=sumaM+eMecProfile(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcretoCirc(a,fdpc,diam);

eleMec=[elemAc;
        elemConc];
      