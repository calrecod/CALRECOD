function eleMec=eleMecanicosRebarTBeams(c,a,fdpc,ha,ba,bp,ht,span,E,rebarType,...
                                    dispositionRebar,rebarAvailable)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicosRebarTBeams(c,a,fdpc,ha,ba,bp,ht,span,E,rebarType,...
%                               dispositionRebar,rebarAvailable)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the sum of resistant forces of a T beam cross-section, 
% considering the contribution of rebars over the cross-section and concrete
% in compression.
% 
% OUTPUT: eleMec:           vector that contains the output [Fs, Ms;
%                                                            Fc, Mc]
%
% INPUT:  rebarType:        Vector that contains the rebar diameters' 
%                           indices of the optimal rebar design (both in 
%                           tension and compression). The vector size is 
%                           nbars x 1 containing a number between 1 
%                           and n-diam according to the available commercial 
%                           rebar types stated by default in the rebar
%                           database
%
%         ha:               flange thickness
%
%         ba:               flange width
%
%         bp:               web width
%
%         ht:               total cross-section height
%
%         c:                neutral axis depth
%
%         a:                reduced neutral axis depth as a = beta1 * c
%
%         dispositionRebar: local coordinates of rebars laid out over the 
%                           beam cross-section
%
%         rebarAvailable:   data base of the available commercial types of 
%                           rebar
%
%         fdpc=0.85f'c:     according to ACI 318
%
%         E:                is the Modulus of Elasticity of the reinforcing
%                           steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-01-19
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

if a<=ha % the beam is working as a T beam
    
    % Effective width of flange be
    be=min([bp+16*ha,bp+ba,bp+span/4]);

    % Centroid of the T cross-section:
    cs=((ha*be)*(0.5*ha)+(ht-ha)*bp*(ha+0.5*(ht-ha)))...
                /(be*ha+(ht-ha)*bp);
            
else % The beam is working as a rectangular beam
    
    % Effective width of flange be
    be=bp;
    
    % Centroid of the effective rectangular cross-section:
    cs=0.5*ht;
    
end

nv=length(rebarType); % total number of rebars

eMecRebar=zeros(nv,9);% [dA, xposition, yposition, d, eps, eE, Fr, Mr]

eMecRebar(:,2)=dispositionRebar(:,1);
eMecRebar(:,3)=dispositionRebar(:,2);


sumaM=0;
sumaF=0;
for i=1:nv
    eMecRebar(i,1)=rebarAvailable(rebarType(i),2)^2*pi/4;
    eMecRebar(i,4)=ht/2-eMecRebar(i,3); % (effective height d)
    eMecRebar(i,5)=0.003/c*(eMecRebar(i,4)-c); % epsilum
    
    if(eMecRebar(i,5)>0.0021)
        eMecRebar(i,5)=0.0021; %compression
    end
    eMecRebar(i,6)=eMecRebar(i,5)*E;
    eMecRebar(i,7)=eMecRebar(i,6)*eMecRebar(i,1);
    eMecRebar(i,8)=eMecRebar(i,7)*(eMecRebar(i,4)-cs);
    
    sumaF=sumaF+eMecRebar(i,7);
    sumaM=sumaM+eMecRebar(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcretoTsec(a,fdpc,bp,be,ha,ht);

eleMec=[elemAc;
        elemConc];
      