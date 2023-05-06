function [raiz]=bisectionMr2tBeamsT(cUno,cDos,fr,E,t1,t2,bp,ht,ba,ha,Lb,...
                                    cover,fdpc,beta1,ea)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMr2tBeamsT(cUno,cDos,fr,E,t1,t2,bp,ht,ba,ha,Lb,cover,...
%                                  fdpc,beta1,ea)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and resistant bending moment 
% of a reinforced T-beam cross-section taking on account both the steel in 
% compression and steel in tension with the aid of the bisection method 
% as a root for the pre-established equilibrium condition sum F=0.
% 
% OUTPUT: raiz:         vector that contains the output [c,sum Fi,MR]
%
% INPUT:  cUno,cDos:    are the initial values for the neutral axis depth 
%                       as the bisection method requires them to begin the 
%                       iterations. Such values are recommended to be 
%                       cUno -> 0 and cDos -> ht
%
%         fr:           is the applied axial force over the beam
%                       cross-section (which for pure flexure is considered
%                       as 0)
%
%         t1,t2         are the given width of ISR in compression and
%                       tension, respectively
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
%         cover:        is the concrete cover for the reinforcing steel
%
%         fdpc:         is the reduced f'c as 0.85f'c according to the
%                       ACI 318-19 code
%
%         beta1:        is determined according to the ACI 318-19 code 
%                       (see documentaiton)
%
%         ea:           is the approximation root error
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-20
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

d1=ht-cover;
%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta1*cUno;
aDos=beta1*cDos;

[eleMec]=eleMecanicos2tBeamsT(cUno,aUno,fdpc,bp,ht,ba,ha,cover,E,t1,t2,Lb);
frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecanicos2tBeamsT(cDos,aDos,fdpc,bp,ht,ba,ha,cover,E,t1,t2,Lb);
frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));

if c==0
    c=0.000001;
elseif c>(d1/(0.005/0.003+1))
    c=(d1/(0.005/0.003+1));
end
a=beta1*c;
[eleMec]=eleMecanicos2tBeamsT(c,a,fdpc,bp,ht,ba,ha,cover,E,t1,t2,Lb);
frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)

    if((raizUno*raizc)<0)
        cDos=c;
        aDos=beta1*cDos;
        [eleMec]=eleMecanicos2tBeamsT(cDos,aDos,fdpc,bp,ht,ba,ha,cover,...
                                      E,t1,t2,Lb);
        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        aUno=beta1*cUno;
        [eleMec]=eleMecanicos2tBeamsT(cUno,aUno,fdpc,bp,ht,ba,ha,cover,...
                                      E,t1,t2,Lb);
        frt=eleMec(1,1)+eleMec(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;

    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c==0
        c=0.000001;
    elseif c>(d1/(0.005/0.003+1))
        c=(d1/(0.005/0.003+1));
    end
    
    a=beta1*c;
    [eleMec]=eleMecanicos2tBeamsT(c,a,fdpc,bp,ht,ba,ha,cover,...
                                  E,t1,t2,Lb);
    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);
    
end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
