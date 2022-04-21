function [maxef,Mrt,c]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
                                     b_rec,h_rec,beta1)

%------------------------------------------------------------------------
% Syntax:
% [maxef,Mrt,c]=Efrec2tBeams(load_conditions,fc,factor_fc,E,h,b,ast,asc,...
%                            b_rec,h_rec,beta1)
%
%------------------------------------------------------------------------
% PURPOSE: Calculates the structural efficiency of beam cross-section 
% according to the applied load conditions.
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
%         h_rec,b_rec:     is the concrete cover along the height dimension
%                          and the width cross-section dimension, respectively
%                          (cm)
%
%         h,b:             cross-section dimensions (cm)
%
%         ast,asc:         are the reinforcement steel area in tension and
%                          compression, respectively (cm2)
%
%         E:               is the Elasticity Modulus of reinforcing steel
%                          (Kg/cm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fdpc=fc*factor_fc;
t1=asc/(b-2*b_rec);
t2=ast/(b-2*b_rec);

d=h-h_rec;
mumax=max(abs(load_conditions(:,2)));

cUno=0.001; 
cDos=h;
fr=0;
[Root]=bisectionMr2tBeams(cUno,cDos,fr,E,t1,t2,h,b,b_rec,h_rec,fdpc,beta1,0.001);  

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


