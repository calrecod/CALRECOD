function [raiz]=bisectionMrRebarTBeams(c1,c2,fr,E,ha,ba,bp,ht,Lb,cover,fdpc,...
    beta1,ea,rebarType,dispositionRebar,rebarAvailable)
%
%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrRebarTBeams(cUno,cDos,fr,E,ha,ba,bp,ht,Lb,cover,fdpc,...
% beta,ea,rebarType,dispositionRebar,rebarAvailable)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and resistant bending moment
% of a reinforced beam T cross-section taking on account the distribution of
% rebars over the cross-section with the aid of the bisection method as a 
% root method for the pre-established equilibrium condition sum F=0.
% 
% OUTPUT: raiz:         vector that contains the neutral axis depth c, the
%                       sum of axial forces of equilibrium sum FR=0 and the
%                       resistant bending moment a [c,sum FR,MR]
%
% INPUT:  ha:           flange thickness
%
%         ba:           flange width
%
%         bp:           web width
%       
%         ht:           total cross-section height
%
%         c1,c2:        initial root values for the use of the bisection 
%                       method. As a closed root method, it is recommended
%                       to use c1=1x10^(-6) and c2=2*ht
%
%         rec:          is the concrete cover parameter vertically, (cm)
%
%         fdpc:         is the reduced f'c as 0.85f'c according to the
%                       ACI 318-19 code
%
%         beta1:        is determined according to the ACI 318-19 code 
%                       (see documentaiton)
%
%         ea:           approximation error to the real value of c
%
%         rebarType:    Vector that contains the rebar diameters' indices 
%                       of the optimal rebar design (both in tension and 
%                       compression) - according to the rebar database 
%                       table. The vector's size is nbars x 1 containing a 
%                       number between 1 and 7 according (by default)
%
%         dispositionRebar: local coordinates of rebars laid out over 
%                           the T-beam cross-section
%
%         rebarAvailable:data base of the available commercial 
%                        rebar sizes. An array of size n# x 2 containing
%                        the eight-of-an-inch available rebars in the first
%                        column and the rebar diameter in the second 
%                        column, from the smallest to the biggest diameter
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-04-08
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

d1=ht-cover;
%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta1*c1;
aDos=beta1*c2;

[eleMec]=eleMecanicosRebarTBeams(c1,aUno,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                                dispositionRebar,rebarAvailable);
frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecanicosRebarTBeams(c2,aDos,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                                dispositionRebar,rebarAvailable);

frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=c2-(raizDos*(c1-c2)/(raizUno-raizDos));

if c==0
    c=0.00001;
elseif c>(d1/(0.005/0.003+1))
    c=(d1/(0.005/0.003+1));
end
a=beta1*c;
[eleMec]=eleMecanicosRebarTBeams(c,a,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                                dispositionRebar,rebarAvailable);

frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=c2;
es=abs((c-cu)/c);
while(es>ea)

    if((raizUno*raizc)<0)
        c2=c;
        aDos=beta1*c2;
        [eleMec]=eleMecanicosRebarTBeams(c2,aDos,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                                    dispositionRebar,rebarAvailable);
        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        c1=c;
        aUno=beta1*c1;
        [eleMec]=eleMecanicosRebarTBeams(c1,aUno,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                                    dispositionRebar,rebarAvailable);
        frt=eleMec(1,1)+eleMec(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;

    c=c2-(raizDos*(c1-c2)/(raizUno-raizDos));
    if c==0
        c=0.000001;
    elseif c>(d1/(0.005/0.003+1))
        c=(d1/(0.005/0.003+1));
    end
    a=beta1*c;
    [eleMec]=eleMecanicosRebarTBeams(c,a,fdpc,ha,ba,bp,ht,Lb,E,rebarType,...
                            	dispositionRebar,rebarAvailable);
    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);

end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
