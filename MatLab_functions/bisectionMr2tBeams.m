function [raiz]=bisectionMr2tBeams(cUno,cDos,fr,E,t1,t2,h,b,b_rec,h_rec,...
                                    fdpc,beta,ea)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMr2tBeams(cUno,cDos,fr,E,t1,t2,h,b,b_rec,h_rec,...
%                           fdpc,beta,ea)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and resistant bending moment 
% of reinforced beam cross-section taking on account both the steel in 
% compression and steel in tension with the aid of the bisection method 
% as a root for the pre-established equilibrium condition sum F=0.
% 
% OUTPUT: raiz:         vector that contains the output [c,sum Fi,MR]
%
% INPUT:  t1,t2         are the given width of ISR in compression and
%                       tension, respectively
%
%         b_rec,h_rec:  are the concrete cover parameters horizontally and
%                       vertically, respectively 
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
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

d1=h-h_rec;
%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta*cUno;
aDos=beta*cDos;

[eleMec]=eleMecanicos2tBeams(cUno,aUno,fdpc,h,b,b_rec,h_rec,E,t1,t2);
frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecanicos2tBeams(cDos,aDos,fdpc,h,b,b_rec,h_rec,E,t1,t2);
frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));

if c==0
    c=0.000001;
elseif c>(d1/(0.005/0.003+1))
    c=(d1/(0.005/0.003+1));
end
a=beta*c;
[eleMec]=eleMecanicos2tBeams(c,a,fdpc,h,b,b_rec,h_rec,E,t1,t2);
frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)

    if((raizUno*raizc)<0)
        cDos=c;
        aDos=beta*cDos;
        [eleMec]=eleMecanicos2tBeams(cDos,aDos,fdpc,h,b,b_rec,h_rec,E,t1,t2);
        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        aUno=beta*cUno;
        [eleMec]=eleMecanicos2tBeams(cUno,aUno,fdpc,h,b,b_rec,h_rec,E,t1,t2);
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
    
    a=beta*c;
    [eleMec]=eleMecanicos2tBeams(c,a,fdpc,h,b,b_rec,h_rec,E,t1,t2);
    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);
    
end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
