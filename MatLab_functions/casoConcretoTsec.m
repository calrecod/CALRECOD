function elemConc=casoConcretoTsec(a,fdpc,bp,be,ha,ht)

%------------------------------------------------------------------------
% Syntax:
% elemConc=casoConcretoTsec(a,fdpc,b,h)
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the contribution of resistance of the concrete
% compression zone of a T reinforced concrete cross-section, regarding 
% axial and bending forces.
% 
% OUTPUT: elemConc: vector that contains the output [F_c, M_c] of 
%                    resistant axial and bending forces
%
% INPUT:  a:        is the reduced depth of neutral axis of the
%                   cross-section in question
%
%         be:       effective flange width in compression
%
%         bp:       web width of the T cross-section
%
%         ha:       flange thickness of the T cross-section
%
%         fdpc:     is the factored value of f'c as 0.85f'c according 
%                   to the ACI 318-19 code  
%
%         ht:       total height of the T beam cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

% Centroid of the T cross-section:
cs=((ha*be)*(0.5*ha)+(ht-ha)*bp*(ha+0.5*(ht-ha)))...
            /(be*ha+(ht-ha)*bp);

if (a>ht)
    a=ht;
end

% area in compression
if a<=ha
    ac=a*be;
    ca=0.5*a;
else
    ac=ha*be+(a-ha)*bp;
    
    % Centroid of the concrete zone in compression
    ca=((ha*be)*(0.5*ha)+(a-ha)*bp*(ha+0.5*(a-ha)))/(be*ha+(a-ha)*bp);
end
frc=-ac*fdpc;
mrc=-frc*(cs-ca); % momentum with respect to cross-section's centroid 

elemConc=[frc mrc];