function elemConc=casoConcretoCirc(a,fdpc,diam)

%------------------------------------------------------------------------
% Syntax:
% elemConc=casoConcretoCirc(a,fdpc,diam)
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the contribution of resistance of the concrete
% compression zone of a circular column cross-section, regarding axial 
% and bending forces.
% 
% OUTPUT: elemConc: vector that contains the output [Fc, Mc] of 
%                    resistant axial and bending forces
%
% INPUT:  a:        is the reduced depth of neutral axis of the
%                   cross-section in question
%
%         fdpc:     is the factored value of f'c as 0.85f'c according 
%                   to the ACI 318-19 code  
%
%         diam:     is the cross-section diameter
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

% computing the compression zone of the concrete
x0=-sqrt((diam*0.5)^2-(a-0.5*diam)^2);
n=100;
dx=(2*abs(x0))/n;
if a<=0.5*diam
    yco=-sqrt((diam*0.5)^2-x0^2)+0.5*diam;
    frc=-fdpc*(a-yco)*dx;
    mrc=fdpc*(a-yco)*dx*(0.5*diam-0.5*(a+yco));
    for i=1:n-1
        x=x0+dx*i;
        y=-sqrt((diam*0.5)^2-x^2)+0.5*diam;
        
        frc=frc-fdpc*(a-y)*dx;
        mrc=mrc+fdpc*(a-y)*dx*(0.5*diam-0.5*(a+y));
    end
    
elseif a>0.5*diam && a<=diam
    yco=sqrt((diam*0.5)^2-x0^2)+0.5*diam;
    frc=-fdpc*(pi*(0.5*diam)^2-(yco-a)*dx);
    
    mrc=fdpc*(-(yco-a)*(0.5*diam-(a+0.5*(yco-a)))*dx);
    for i=1:n-1
        x=x0+dx*i;
        y=sqrt((diam*0.5)^2-x^2)+0.5*diam;
        
        frc=frc+fdpc*(y-a)*dx;
        mrc=mrc+fdpc*(-(y-a)*dx*(0.5*diam-(a+0.5*(y-a))));
    end
elseif (a>diam)
    frc=fdpc*(-pi*(0.5*diam)^2);
    mrc=0;
end

elemConc=[frc mrc];

