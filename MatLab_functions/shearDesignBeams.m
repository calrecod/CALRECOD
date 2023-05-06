function [s1,s2,s3,d1,d2]=shearDesignBeams(span,b,h,rec,rho,fc,fy,...
                          shear_beam)
                      
%------------------------------------------------------------------------
% Syntax:
% [s1,s2,s3,d1,d2]=shearDesignBeams(span,b,h,rec,rho,fc,fy,shear_beam)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%------------------------------------------------------------------------
% PURPOSE: To design the separation of the transversal reinforcement along
% the whole length of a beam element according to the mechanic shear
% forces distribution.
% 
% OUTPUT: s1:                   is the separation of the transversal 
%                               reinforcement at the left part/section of
%                               the beam element.
%
%         s2:                   is the separation of the transversal 
%                               reinforcement at the middle part/section 
%                               of the beam element.
%
%         s3:                   is the separation of the transversal 
%                               reinforcement at the right part/section 
%                               of the beam element.
%                               
%         d1:                   is the length along which the transversal
%                               reinforcement is separated by s1 (cm).
%                               Left part length.
%                               
%         d2:                   is the length along which the transversal
%                               reinforcement is separated by s3 (cm).
%                               Right part length.
%
% INPUT:  span:                 is the length of the beam element
%
%         b,h:                  are the cross-section dimensions: width 
%                               and height, respectively.
%
%         rec:                  is the concrete cover
%
%         rho:                  is the average longitudinal reinforcing 
%                               percentage area
%
%         shear_beam:           is the array containing the shear forces
%                               distribution from left to right.
%
%         fc, fy:               is the compressive strength of the concrete
%                               and the yield stress of the reinforcing
%                               steel.
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-03-19
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

np=length(shear_beam);

%% Cross-section shear resistance of the beam (Vcr)
Frv=0.75; % shear resistance factor
Vcr=Frv*(b*(h-rec))*(0.2+20*rho)*sqrt(fc);

%% Stirrups separation
diam=3/8*2.54;
ab=diam^2*pi*0.25;
% Left part
s1=abs((Frv*2*ab*fy*(h-rec))/(abs(shear_beam(1))-Vcr));
if s1<6
    s1=6;
elseif s1>(h-rec)*0.5
    s1=(h-rec)*0.5;
end
s1=s1-mod(s1,5);

% Right part
s3=abs((Frv*2*ab*fy*(h-rec))/(abs(shear_beam(np))-Vcr));
if s3<6
    s3=6;
elseif s3>(h-rec)*0.5
    s3=(h-rec)*0.5;
end
s3=s3-mod(s3,5);

% Middle part
s2=(s1+s3)*0.5;
s2=s2-mod(s2,5);
Vus2=Vcr+abs((Frv*2*ab*fy*(h-rec))/s2);

% To determine the position in which the separation of the stirrups 
% will change from s1 to s2 along the beam beam length, starting from
% the left
ps21=[];
for i=1:np-1
    if abs(shear_beam(i))<=Vus2 && abs(shear_beam(i+1))>=Vus2
        ps21=i;
    end
end

ps22=[];
for i=1:np-1
    if abs(shear_beam(np+1-i))>=Vus2 && abs(shear_beam(np-i))<=Vus2
        ps22=np-i;
    end
end

lenbeam=linspace(0,span,np);
if isempty(ps21)==0
    d1=lenbeam(ps21);
else
    d1=ceil(span/2);
   
end
if isempty(ps22)==0
    d2=span-lenbeam(ps22);
else
    d2=span-d1;
end