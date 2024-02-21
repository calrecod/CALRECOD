function [root]=bisectionMrSymISRRot(cUno,cDos,fr,E,h,b,fdpc,beta1,...
        ea,ndt,da4,coordiscISR,cp,RotCornerSec,gamma)
%------------------------------------------------------------------------
% Syntax:
% [root]=bisectionMrSymISRRot(cUno,cDos,fr,E,h,b,fdpc,beta1,...
%      ea,ndt,da4,coordiscISR,cp,RotCornerSec,gamma)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and bending resistance
% from the interaction diagram of a rotated reinforced concrete column 
% cross-section reinforced with a symmetrical ISR, corresponding to a
% certain axial load resistance with the aid of the bisection root method.
% 
% OUTPUT: root:                 is a vector containing the neutral axis 
%                               depth, axial resistant force and bending 
%                               resistance of a reinforced column cross
%                               section as [c,FR,MR] (neutral axis depth, 
%                               resisting axial forces and resisting moment)
%
% INPUT:  cUno,cDos:            are the initial values of the neutral axis
%                               to commence iterations
%
%         fr:                   is the axial force resistance corresponding
%                               to the bending moment resistance for which
%                               the equilibrium condition sum F=0 is 
%                               established to extract its corresponding 
%                               neutral axis depth from the interaction
%                               diagram. For columns, the value of this 
%                               variable may go from the max resistance 
%                               in compression of the cross-section (poc)
%                               to the max resistance in tension (pot)
%
%         E:                    Elasticity modulus of steel
%
%         b,h:                  cross-section dimensions
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta1:                 is determined as stablished by code 
%                               (see Documentation)
%
%         ea:                   is the approximation error to terminate the
%                               root bisection method
%
%         ndt:                  is the total number of ISR's discrete
%                               elements over the cross-section
%
%         coordiscISR:          are the local coordinates of the discrete
%                               ISR's elements over the cross-section
%
%         cp:                   Plastic Center location for each axis
%                               direction of the rotated column 
%                               cross-section
%
%         RotCornerSec:         are the coordinates of each of the four
%                               corners of the rotated cross-section
%
%         gamma:                is the angle of rotation of the
%                               column cross-section
%
%         da4:                  are the ISR's discrete elements'
%                               cross-section area, for each of the four
%                               column cross-section's boundaries
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%

[eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,cUno,fdpc,...
    E,beta1,cp,RotCornerSec,gamma);
                              
frt=eMecISR(1,1)+eMecISR(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,cDos,fdpc,...
    E,beta1,cp,RotCornerSec,gamma);

frt=eMecISR(1,1)+eMecISR(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));

if c==0
    c=1e-6;
end

[eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,c,fdpc,...
    E,beta1,cp,RotCornerSec,gamma);

frt=eMecISR(1,1)+eMecISR(2,1);
raizc=fr-frt;

%% Main loop
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)
    
    if((raizUno*raizc)<0)
        cDos=c;
        [eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,cDos,...
            fdpc,E,beta1,cp,RotCornerSec,gamma);

        frt=eMecISR(1,1)+eMecISR(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        [eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,cUno,...
            fdpc,E,beta1,cp,RotCornerSec,gamma);

        frt=eMecISR(1,1)+eMecISR(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;

    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c==0
        c=0.00001;
    end
    [eMecISR]=eleMecISRSymRecRot(coordiscISR,ndt,da4,b,h,c,fdpc,E,...
        beta1,cp,RotCornerSec,gamma);

    frt=eMecISR(1,1)+eMecISR(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);

    if itdos>100 || ituno>100
        break;
    end
end
mrt=eMecISR(1,2)+eMecISR(2,2);
root=[c,frt,mrt];   