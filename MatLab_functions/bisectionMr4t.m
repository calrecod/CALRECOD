function [raiz]=bisectionMr4t(cUno,cDos,fr,E,t1,t2,t3,t4,h,b,rec,fdpc,beta1,ea)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMr4t(cUno,cDos,fr,E,t1,t2,t3,t4,h,b,rec,fdpc,beta1,ea)
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and bending moment resistance
% from the interaction diagram of a reinforced column cross-section given 
% a resistant axial load.
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
%         t1,t2,t3,t4:        are the ISR depths of a 4t-ISR (strictly for
%                             columns)
%
%         rec:                is a vector containing the concrete cover for
%                             both cross-section axis as [coverX,coverY]
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code
%
%         beta1:              is determined as established in ACI 318 code
%                             (see Documentation)
%
%         ea:                 is the approximation root error (close to 0)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta1*cUno;
aDos=beta1*cDos;

[eleMec]=eleMecanicos4t(cUno,aUno,fdpc,h,b,rec,E,t1,t2,t3,t4);
frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecanicos4t(cDos,aDos,fdpc,h,b,rec,E,t1,t2,t3,t4);
frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
if c==0
    c=0.00000001;
end
a=beta1*c;

[eleMec]=eleMecanicos4t(c,a,fdpc,h,b,rec,E,t1,t2,t3,t4);
frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%% inicia ciclo %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
total_iterations=0;
while(es>ea)

    if((raizUno*raizc)<0)
        cDos=c;
        aDos=beta1*cDos;
        [eleMec]=eleMecanicos4t(cDos,aDos,fdpc,h,b,rec,E,t1,t2,t3,t4);
        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        aUno=beta1*cUno;
        [eleMec]=eleMecanicos4t(cUno,aUno,fdpc,h,b,rec,E,t1,t2,t3,t4);
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
    [eleMec]=eleMecanicos4t(c,a,fdpc,h,b,rec,E,t1,t2,t3,t4);
    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);
    
    if total_iterations>1500
        break;
    end
    
    total_iterations=total_iterations+1;
    
end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
