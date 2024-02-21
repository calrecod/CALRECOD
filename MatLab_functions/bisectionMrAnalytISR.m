function [raiz]=bisectionMrAnalytISR(cUno,cDos,fr,E,t,h,b,rec,fdpc,beta1,ea)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrAnalytISR(cUno,cDos,fr,E,t,h,b,rec,fdpc,beta1,ea)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and bending moment resistance
% from the interaction diagram of a rectangular reinforced column 
% cross-section given a corresponding resistant axial load. The bisection
% method is used.
% 
% OUTPUT: raiz:               vector containing the neutral axis depth, 
%                             axial resistant force and bending resistance 
%                             of a reinforced column cross-section as 
%                             [c,FR,MR]
%
% INPUT:  cUno,cDos:          Are the initial neutral axis depth values for
%                             the implementation of the bisection method,
%                             recommended to be close to the values:
%                             cUno=0.0001 and cDos=2h
%
%         fr:                 is the axial load resistance from which the
%                             bending resistance will be determined from 
%                             the interaction diagram
%
%         t:                  is the ISR width of a 1t-ISR (strictly for
%                             columns)
%
%         rec:                is the concrete cover
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code
%                             fdpc=0.85*fc
%
%         beta1:              is determined as established in ACI 318 code
%                             (see Documentation)
%
%         ea:                 is the approximation root error (close to 0)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta1*cUno;
aDos=beta1*cDos;

rec=rec(1);
dUno=rec+0.5*t;
dDos=h-rec-0.5*t;
[eleMec]=eleMecISRAnalyt(cUno,aUno,fdpc,h,b,rec,E,t,dUno,dDos);

frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecISRAnalyt(cDos,aDos,fdpc,h,b,rec,E,t,dUno,dDos);

frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
if c==0
    c=0.00000001;
end
a=beta1*c;
[eleMec]=eleMecISRAnalyt(c,a,fdpc,h,b,rec,E,t,dUno,dDos);

frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
total_iterations=0;
while(es>ea)

    if((raizUno*raizc)<0)
        cDos=c;
        aDos=beta1*cDos;
        [eleMec]=eleMecISRAnalyt(cDos,aDos,fdpc,h,b,rec,E,t,dUno,dDos);

        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        aUno=beta1*cUno;
        [eleMec]=eleMecISRAnalyt(cUno,aUno,fdpc,h,b,rec,E,t,dUno,dDos);

        frt=eleMec(1,1)+eleMec(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;
    
    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c<=0
        c=0.001;
    end
        
    a=beta1*c;
    [eleMec]=eleMecISRAnalyt(c,a,fdpc,h,b,rec,E,t,dUno,dDos);

    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);
    
    if total_iterations>100
        break;
    end
    
    total_iterations=total_iterations+1;
    
end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
