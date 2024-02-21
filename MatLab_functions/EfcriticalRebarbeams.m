function [maxef,Mrv,c]=EfcriticalRebarbeams(load_conditions,b,E,fdpc,arrange_t1,...
    arrange_t2,rebarAvailable,d,h_rec,beta1,disposition_rebar)

%------------------------------------------------------------------------
% Syntax:
% [maxef,Mrv,c]=EfcriticalRebarbeams(load_conditions,b,E,fdpc,arrange_t1,...
%   arrange_t2,rebarAvailable,d,h_rec,beta1,disposition_rebar)
%
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the bending resistance, structural efficiency and 
% depth of neutral axis for the optimal reinforced beam cross-section with
% rebars.
% 
% OUTPUT: maxef:                structural efficiency for the optimal 
%                               reinforced cross-section with rebars
%
%         Mrv:                  Resistant bending moment for the optimal 
%                               reinforced cross-section with rebars
%
%         c:                    neutral axis depth for the optimal 
%                               reinforced cross-section with rebars
%
% INPUT:  h_rec:                is the concrete cover vertically
%
%         fdpc:                 is the f'c reduced by the factor 0.85 
%                               according to code
%
%         E:                    is the Elasticity Modulus of steel
%
%         load_conditions:      Is the array containing the load conditions
%                               in format: [nload, Mu]
%
%         rebarAvailable:       Commercial available type of rebars database
%                               as indicated in function:
%                               "ISR1tRebarBeamsOptimization"
%                               (see Documentation)
%
%         disposition_rebar:    local coordinates of rebars laid out over 
%                               the beam cross-section
%
%         arrange_t1,           Vectors that contain the type of rebar for
%         arrange_t2:           the optimal option both in tension and
%                               compression, respectively. The vectors size
%                               is of one column with nrebar rows containing
%                               a number between 1 and 7 according to the 
%                               available commercial rebar types stated by 
%                               default
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

h=d+h_rec;
mumax=max(abs(load_conditions(:,2)));

cUno=0.001; 
cDos=h;
fr=0;

[Root]=bisectionMrRebarBeams(cUno,cDos,fr,E,h,b,h_rec,fdpc,beta1,0.001,...
    arrange_t1,arrange_t2,disposition_rebar,rebarAvailable);  

c=Root(1);
Mrv=Root(3);

eps_tens_ac=0.003/c*(d-c);

if eps_tens_ac>0.005
    factor_resistance=0.9;
else
    factor_resistance=0.65+(eps_tens_ac-0.0021)*250/3;
end
Mrv=Mrv*factor_resistance;
maxef=mumax/Mrv;