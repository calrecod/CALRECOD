function [maxef,mr]=EvaluateISR1tFoot(t_tension,b,h,...
                fy,fdpc,rec,beta1,axis,mu_real_axis)

%------------------------------------------------------------------------
% Syntax:
% [maxef,mr]=EvaluateISR1tFoot(t_tension,b,h,...
%                fy,fdpc,rec,beta1,axis,mu_real_axis)
%
%------------------------------------------------------------------------
% PURPOSE: To determine the structural efficiency of a footing transversal
% cross-section subject to pure flexure.
% 
% OUTPUT: maxef:                is the structural efficiency of the 
%                               reinforced footing transversal cross-section
%                               as maxef=Mu/MR
%
%         mr:                   is the resistant bending moment of the 
%                               reinforced footing transversal cross-section
%                               MR
%
% INPUT:  h:                    is the cross-section height
%
%         b:                    is the cross-section width dimension
%
%         fdpc:                 is the reduced f'c as fdcp=0.85f'c according
%                               to code
%
%         rec                   is the concrete cover
%
%         mu_real_axis:         is the effective flexure distributed to the
%                               transversal cross-section from the contact
%                               soil pressures
%
%         fy:                   is the yield stress of reinforcing steel
%
%         axis:                 is the footing axis direction of analysis:
%                               (1) represents the axis direction in which 
%                               the dimension L is the width of the 
%                               transversal cross-section, for (2) B is the
%                               width of the transversal cross-section, in
%                               its reference system (see Documentation)
%
%         beta1:                is determined as specified in code (see 
%                               Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bp=b-2*rec;

ac_tension=t_tension*bp;
            
if axis==1
    d=h-rec;
elseif axis==2
    d=h-rec-12/8*2.54;
end
    
c=ac_tension*fy/(beta1*b*fdpc);
if c>(d/(0.004/0.003+1))
    c=d/(0.004/0.003+1);
end
a=0.85*c;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Efficiency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=ac_tension*fy;

mr=0.9*(T*(d-a/2))*1e-5;
maxef=mu_real_axis/mr;

       