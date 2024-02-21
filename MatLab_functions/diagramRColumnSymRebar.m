function [diagrama,cPoints,poc,pot]=diagramRColumnSymRebar(As,b,h,E,npuntos,...
                       fdpc,nv,beta,ov,av,rebarDisposition)
                            
%------------------------------------------------------------------------
% Syntax:
% [diagrama,cPoints,poc,pot]=diagramRColumnSymRebar(As,b,h,E,npuntos,...
% fdpc,nv,beta,ov,av,disposicion_varillado)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram of a symmetric rebar option, 
% as well as the structural efficiency given certain load conditions.
% 
% OUTPUT: diagrama:             is the interaction diagram data for both
%                               cross-section's axis. Format: 
%
%               [P,MRx,FR*P,FR*MRx,ec-x,MRy,FR*P,FR*MRy,ecc-y]
%
%         cPoints:              are the neutral axis depth values 
%                               corresponding to each point of the 
%                               interaction diagram
%
% INPUT:  b,h:                  given cross-section dimensions
%
%         av,ov,nv:             are the rebar area, the type of rebar in 
%                               eighth of inches (ov/8 in) and the number
%                               of rebars, respectively
%
%         beta:                 is determined as specified in code ACI,
%                               according to the f'c used
%
%         fdpc:                 equal to 0.85f'c according to code
%
%         npuntos:              number of points to compute for the
%                               definition of the interaction diagram
%
%         rebarDisposition:     Local coordinates of the rebars over the
%                               cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

if npuntos<3
    disp('Error: the number of points for the Interaction Diagram must be 3 or higher');
    return;
end
    
diagrama=zeros(npuntos,9);
cPoints=zeros(npuntos,2);

fy=4200;
Ac=b*h;
pot=As*fy;

poc=(fdpc*(Ac-As)+fy*As);
df=(pot+poc)/(npuntos-1);

diagrama(1,1)=-poc;
diagrama(npuntos,1)=pot;
diagrama(1,2)=0;
diagrama(1,6)=0;

cPoints(1,1)=4*h;
cPoints(1,2)=4*h;
dfi=0.1/(npuntos-1);

varillado=[rebarDisposition(:,1) rebarDisposition(:,2)];
dimensionesColumna=[b h];
for sentido=1:2
    if (sentido==2)
        b=dimensionesColumna(2);
        h=dimensionesColumna(1);
        rebarDisposition(:,1)=-varillado(:,2);
        rebarDisposition(:,2)=varillado(:,1);
    end
    diagrama(1,4*sentido-1)=0.65*-poc;
    diagrama(1,4*sentido)=0;
    for i=1:npuntos-1
        diagrama(i+1,1)=(-poc+i*df);
        fr=diagrama(i+1,1);
        cUno=0.001;
        cDos=4*h;

        [raiz]=bisectionMrSymRebarCols(cUno,cDos,fr,E,h,b,fdpc,beta,...
                                 0.001,nv,ov,av,rebarDisposition);  

        diagrama(i+1,4*sentido-2)=raiz(3);

        cPoints(i+1,sentido)=raiz(1);

        %%%%%%%%%%%%%%% Reduced diagrams %%%%%%%%%%%%%%%%%%%%
        diagrama(i+1,4*sentido-1)=(0.65+i*dfi)*diagrama(i+1,1);
        diagrama(i+1,4*sentido)=(0.65+i*dfi)*diagrama(i+1,4*sentido-2);
        
        %%%%%%%%%%%%%%%%%%%%%%%% Eccentricities %%%%%%%%%%%%%%%%%%%%%%%

       diagrama(i+1,4*sentido+1)=diagrama(i+1,4*sentido)/...
                                         diagrama(i+1,4*sentido-1);
    end
    
    diagrama(npuntos,4*sentido-1)=0.75*pot;
    diagrama(npuntos,4*sentido)=0;
    diagrama(npuntos,4*sentido-2)=0;
end
