function elemConc=casoConcreto(a,fdpc,b,h)

%------------------------------------------------------------------------
% Syntax:
% elemConc=casoConcreto(a,fdpc,b,h)
%------------------------------------------------------------------------
% PURPOSE: To compute the contribution of resistance of the concrete
% compression zone of a rectangular beam cross-section, regarding axial 
% and bending forces.
% 
% OUTPUT: elemConc: vector that contains the output $F_c, M_c]$ of 
%                    resistant axial and bending forces
%
% INPUT:  a:        is the reduced depth of neutral axis of the
%                   cross-section in question
%
%         fdpc:     is the factored value of f'c as 0.85f'c according 
%                   to the ACI 318-19 code  
%
%         b,h:      are the cross-section dimensions
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

if (a>h)
    a=h;
end
frc=-a*b*fdpc*0.001;
mrc=-frc*(0.5*h-a/3)*0.01;

elemConc=[frc mrc];
