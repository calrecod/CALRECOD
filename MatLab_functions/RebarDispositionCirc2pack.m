function [dispositionRebar,separation]=RebarDispositionCirc2pack(diam,rec,dv,nv)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar,separation]=RebarDispositionCirc2pack(diam,rec,dv,nv)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the local position coordinates of a rebar design,
%          symmetrically over a circular column cross-section, in packages
%          of two.
% 
% OUTPUT: dispositionRebar:     are the local position coordinates of the 
%                               symmetrical rebar option
%
% INPUT:  diam:                 cross-section diameter
%
%         dv,nv:                are the rebar diameter of the rebar
%                               option and the number of rebars,
%                               respectively
%
%         rec:                  is the concrete cover
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

separation=round((pi*(diam-2*rec)-nv*dv)/nv,1);
nv2=nv*0.5;
nv4=0.5*nv2;
%% Outer rebar
dispositionRebar=zeros(nv,2);
teta=360/nv2;
for j=1:nv2
    x=(diam-2*rec)*0.5*cos(deg2rad(teta*0.5+teta*(j-1)));
    dispositionRebar(j,1)=x;
    
    if j<=nv4
        y=sqrt(((diam-2*rec)*0.5)^2-x^2);
        dispositionRebar(j,2)=y;
    elseif j>nv4
        y=-sqrt(((diam-2*rec)*0.5)^2-x^2);
        dispositionRebar(j,2)=y;
    end

end

%% Inner rebar

for j=1:nv2
    x=(diam-2*rec-2*dv)*0.5*cos(deg2rad(teta*0.5+teta*(j-1)));
    dispositionRebar(nv2+j,1)=x;
    
    if j<=nv4
        y=sqrt(((diam-2*rec-2*dv)*0.5)^2-x^2);
        dispositionRebar(nv2+j,2)=y;
    elseif j>nv4
        y=-sqrt(((diam-2*rec-2*dv)*0.5)^2-x^2);
        dispositionRebar(nv2+j,2)=y;
    end

end
