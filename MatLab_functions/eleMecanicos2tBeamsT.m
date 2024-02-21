function eleMec=eleMecanicos2tBeamsT(c,a,fdpc,bp,ht,ba,ha,cover,Es,...
                                    t1,t2,span)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicos2tBeamsT(c,a,fdpc,bp,ht,ba,ha,cover,Es,t1,t2)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the sum of resistant forces of a T-beam 
% cross-section considering the contribution of steel in tension, steel in 
% compression and concrete zone in compression.
% 
% OUTPUT: eleMec:       vector that contains the output [Fs, Ms;
%                                                        Fc, Mc]
%
% INPUT:  t1,t2         are the given width of ISR in compression and 
%                       tension, respectively
%
%         cover:        is the concrete cover horizontally and vertically
%
%         fdpc:         Is the reduced f'c as 0.85f'c according to ACI 318
%
%         c:            is the nuetral axis depth value
%
%         a:            is the reduced effective neutral axis depth value
%                       a = beta1 * c
%
%         ba:           is the effective flange width of the T-beam 
%                       cross-section
%
%         ht:           is total height of the T-beam cross-section
%
%         bp:           is the web width of the T-beam cross-section
%
%         ha:           is the flange thickness of the T-beam
%                       cross-section
%
%         Lb:           is the length of the beam element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

bpp=bp-2*cover;

eMecProfile=zeros(2,8);

% eMecProfile=[dA, xposition, yposition, d, eps, eE, Fr, Mr]

%calculate dA for each t area
A1=bpp*t1;
A2=bpp*t2;

% localtion of each discrete dt over the section
for i=1:2
    if i==1 % compression ISR
        eMecProfile(i,1)=A1;
        x=0; eMecProfile(i,2)=x;
        
        y=0.5*ht-cover-0.5*t1;
        eMecProfile(i,3)=y;
        
    elseif i==2 % tension ISR
        eMecProfile(i,1)=A2;
        x=0;
        eMecProfile(i,2)=x;
        
        y=-0.5*ht+cover+0.5*t2;
        eMecProfile(i,3)=y;
    end
end

if a<=ha % the beam is working as a T beam
    
    % Effective width of flange be
    be=min([bp+16*ha,bp+ba,bp+span/4]);

    % Centroid of the T cross-section:
    cs=((ha*be)*(0.5*ha)+(ht-ha)*bp*(ha+0.5*(ht-ha)))...
                /(be*ha+(ht-ha)*bp);
            
else % The beam is working as a rectangular beam
    
    % Effective width of flange (be) equal to the web's width
    be=bp;
    
    % Centroid of the effective rectangular cross-section:
    cs=0.5*ht;
    
end

sumaM=0;
sumaF=0;
for i=1:2
    eMecProfile(i,4)=0.5*ht-eMecProfile(i,3); %d1
    eMecProfile(i,5)=0.003/c*(eMecProfile(i,4)-c); %epsilum
    
    if(eMecProfile(i,5)>0.0021)
        eMecProfile(i,5)=0.0021; % Only the steel in compression
                                 % is limited by its yielding strain
    end
    eMecProfile(i,6)=eMecProfile(i,5)*Es;
    eMecProfile(i,7)=eMecProfile(i,6)*eMecProfile(i,1);
    eMecProfile(i,8)=eMecProfile(i,7)*(eMecProfile(i,4)-cs);
    
    sumaF=sumaF+eMecProfile(i,7);
    sumaM=sumaM+eMecProfile(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcretoTsec(a,fdpc,bp,be,ha,ht);

eleMec=[elemAc;
        elemConc];
