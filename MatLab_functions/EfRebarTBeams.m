function [Eff,Mr,c]=EfRebarTBeams(load_conditions,bp,ht,ba,ha,Lb,fdpc,...
                                   rebarType,rebarAvailable,cover,beta1,...
                                   rebarDisposition)
%------------------------------------------------------------------------
% Syntax:
% [Eff,Mr,c]=EfRebarTBeams(load_conditions,bp,ht,ba,ha,Lb,fdpc,...
% rebarType,rebarAvailable,cover,beta1,rebarDisposition)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: Computes the structural efficiency of a rebar reinforced 
% T-beam cross-section according to the applied load conditions (pure 
% flexure).
% 
% OUTPUT: Eff,Mr,c:        is the structural efficiency of the reinforced 
%                          beam cross-section (0-1), the resistant beding
%                          moment and the neutral axis depth, respectively
%
% INPUT:  load_conditions: vector as [nload,Mu] size: nloads x 2 
%
%         beta1:           is determined as preescribed by the ACI 318 code
%                          (see documentation) according to the f'c value
%
%         fdpc:            is the reduced concrete compressive strength 
%                          (0.85 * f'c ) as prescribed in the ACI 318 code
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
%         cover:           is the concrete cover for the reinforcement
%
%         rebarAvailable:  is the rebar data base table with each of the
%                          eight-of-an-inch rebars available (from the
%                          smallest to the biggest diameter)
%
%         rebarType:       is the vector containing the rebar diameters'
%                          indices (according to their place in the rebar 
%                          database table)
%
%         rebarDisposition:is the array containing the rebar local
%                          coordinates over the T-beam cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-04-08
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

if fc<2000 % units: (kg,cm)
    Es=2.0e6; % Modulus of elasticity of rebar (Kg/cm2)
else       % units: (lb,in)
    Es=60000/0.0021; % Modulus of elasticity of rebar (psi)
end
d=ht-cover;
Mu=max(load_conditions(:,2));

cUno=0.001; 
cDos=ht;
fr=0;
[Root]=bisectionMrRebarTBeams(cUno,cDos,fr,Es,ha,ba,bp,ht,Lb,cover,fdpc,...
    beta1,0.001,rebarType,rebarDisposition,rebarAvailable);  

c=Root(1);
Mrv=Root(3);

eps_tens_ac=0.003/c*(d-c);

if eps_tens_ac>0.005
    factor_resistance=0.9;
else
    factor_resistance=0.65+(eps_tens_ac-0.0021)*250/3;
end

Mr=Mrv*factor_resistance;
Eff=Mu/Mr;