function Inertia_modif=InertiaBeamCrackedSection(fc,E,areabartension,...
                b,h,h_rec)
%------------------------------------------------------------------------
% Syntax:
% Inertia_modif=InertiaBeamCrackedSection(fc,E,area_var_tension,...
%               b,h,h_rec)
%------------------------------------------------------------------------
% PURPOSE: To determine the modified inertia momentum of a cracked beam 
%          cross-section
% 
% OUTPUT: Inertia_modif:  The modified inertia momentum for the beam 
%                         cross-section, considering it is cracked, that 
%                         is, when the tension stress is greater than the 
%                         rupture modulus of the used concrete.
%
% INPUT:  areabartension: rebar area in tension
%
%         fc:             Is the f'c used
%
%         E:              Is the Elasticity Modulus (Kg/cm2)  
%
%         b,h:            are the cross-section dimensions
%         
%         h_rec           is the concrete cover along the height dimension
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

d=h-h_rec;

%%% modified intertia_________________
n=E/(14000*(fc)^0.5);

eje_neutro=(-n*areabartension(1)+sqrt((n*areabartension(1))^2+2*b*n*...
            areabartension(1)*d))/b;
I_ag=b*eje_neutro^3/12+b*eje_neutro^3/4+n*areabartension(1)*(d-eje_neutro)^2;
Inertia_modif=I_ag;