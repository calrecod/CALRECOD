function [raiz]=bisectionMrRebarBeams(c1,c2,fr,E,h,b,h_rec,fdpc,beta,...
            ea,arreglo_t1,arreglo_t2,dispositionRebar,rebarAvailable)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrRebarBeams(cUno,cDos,fr,E,h,b,h_rec,fdpc,beta,...
%            ea,arreglo_t1,arreglo_t2,dispositionRebar,rebarAvailable)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth and resistant bending moment
% of a reinforced beam cross-section taking on account the distribution of
% rebars over the cross-section with the aid of the bisection method as a 
% root for the pre-established equilibrium condition sum F=0.
% 
% OUTPUT: raiz:         vector that contains the neutral axis depth c, the
%                       sum of axial forces of equilibrium sum FR=0 and the
%                       resistant bending moment a [c,sum FR,MR]
%
% INPUT:  c1,c2:        initial root values for the use of the bisection 
%                       method. As a closed root method, it is recommended
%                       to use c1=1x10^(-6) and c2=2h
%
%         h_rec:        is the concrete cover parameter vertically
%
%         fdpc:         is the reduced f'c as 0.85f'c according to the
%                       ACI 318-19 code
%
%         beta:         is determined according to the ACI 318-19 code 
%                       (see documentaiton)
%
%         arreglo_t1,   Vectors that contain the type of rebar for the 
%         arreglo_t2:   optimal option both in tension and compression,
%                       respectively. The vectors size is of one column 
%                       with nrebar rows containing a number between 1 and 
%                       n# according to the stated available commercial  
%                       rebar diameters
%
%         dispositionRebar:      local coordinates of rebars laid out over 
%                                the beam cross-section
%
%         rebarAvailable:        data base of the available commercial 
%                                types of rebar
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

d1=h-h_rec;
%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%
aUno=beta*c1;
aDos=beta*c2;

[eleMec]=eleMecanicosRebarBeams(c1,aUno,fdpc,h,b,h_rec,E,arreglo_t1,...
                        arreglo_t2,dispositionRebar,rebarAvailable);
frt=eleMec(1,1)+eleMec(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eleMec]=eleMecanicosRebarBeams(c2,aDos,fdpc,h,b,h_rec,E,arreglo_t1,...
                        arreglo_t2,dispositionRebar,rebarAvailable);

frt=eleMec(1,1)+eleMec(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=c2-(raizDos*(c1-c2)/(raizUno-raizDos));

if c==0
    c=0.00001;
elseif c>(d1/(0.005/0.003+1))
    c=(d1/(0.005/0.003+1));
end
a=beta*c;
[eleMec]=eleMecanicosRebarBeams(c,a,fdpc,h,b,h_rec,E,arreglo_t1,...
                        arreglo_t2,dispositionRebar,rebarAvailable);

frt=eleMec(1,1)+eleMec(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=c2;
es=abs((c-cu)/c);
while(es>ea)

    if((raizUno*raizc)<0)
        c2=c;
        aDos=beta*c2;
        [eleMec]=eleMecanicosRebarBeams(c2,aDos,fdpc,h,b,h_rec,E,arreglo_t1,...
                            arreglo_t2,dispositionRebar,rebarAvailable);
        frt=eleMec(1,1)+eleMec(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        c1=c;
        aUno=beta*c1;
        [eleMec]=eleMecanicosRebarBeams(c1,aUno,fdpc,h,b,h_rec,E,arreglo_t1,...
                            arreglo_t2,dispositionRebar,rebarAvailable);
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
    a=beta*c;
    [eleMec]=eleMecanicosRebarBeams(c,a,fdpc,h,b,h_rec,E,arreglo_t1,...
                            arreglo_t2,dispositionRebar,rebarAvailable);
    frt=eleMec(1,1)+eleMec(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);
    if itdos>100
        break;
    elseif ituno>100
        break;
    end
end
mrt=eleMec(1,2)+eleMec(2,2);
raiz=[c,frt,mrt];   
