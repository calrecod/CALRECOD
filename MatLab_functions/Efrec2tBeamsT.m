function [maxef,Mrt,c]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
                        ha,Lb,ast,asc,cover,beta1)

%------------------------------------------------------------------------
% Syntax:
% [maxef,Mrt,c]=Efrec2tBeamsT(load_conditions,fc,factor_fc,E,bp,ht,ba,...
%                       ha,Lb,ast,asc,cover,beta1)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: Computes the structural efficiency of a T-beam cross-section 
% according to the applied load conditions (pure flexure).
% 
% OUTPUT: maxEf,Mr,c:      is the structural efficiency of the reinforced 
%                          beam cross-section (0-1), the resistant beding
%                          moment and the neutral axis depth, respectively
%
% INPUT:  load_conditions: vector as [nload,Mu] size: nloads x 2
%
%         factor_fc:       is determined by de applicable design code. 
%                          The ACI 318-19 specifies it as 0.85
%
%         beta1:           is determined as preescribed by the ACI 318 code
%                          (see documentation)
%
%         cover:           is the concrete cover for the reinforcement
%
%         ba:              is the effective flange width of the T-beam 
%                          cross-section
%
%         ht:              is total height of the T-beam cross-section
%
%         bp:              is the web width of the T-beam cross-section
%
%         ha:              is the flange thickness of the T-beam
%                          cross-section
%
%         Lb:              is the length of the beam element
%
%         ast,asc:         are the reinforcement steel area in tension and
%                          compression, respectively 
%
%         E:               is the Elasticity Modulus of reinforcing steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

fdpc=fc*factor_fc;
t1=asc/(bp-2*cover);
t2=ast/(bp-2*cover);

d=ht-cover;
mumax=max(abs(load_conditions(:,2)));

cUno=0.001; 
cDos=ht;
fr=0;
[Root]=bisectionMr2tBeamsT(cUno,cDos,fr,E,t1,t2,bp,ht,ba,ha,Lb,cover,...
       fdpc,beta1,0.001);  

c=Root(1);
Mrt=Root(3);

eps_tens_ac=0.003/c*(d-c);
if eps_tens_ac>=0.005
    factor_resistance=0.9;
else
    factor_resistance=0.65+(eps_tens_ac-0.0021)*250/3;
end
Mrt=Mrt*factor_resistance;

maxef=mumax/Mrt;
% ----------------------------------------------------------------------