function [root]=bisectionMrRebarColsCir(cUno,cDos,fr,E,diam,fdpc,beta1,...
                             ea,av,rebarDisposition,rec)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrRebarColsCir(cUno,cDos,fr,E,diam,fdpc,beta1,...
%           ea,av,rebarDisposition,rec)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth, axial and bending resistance
% from the interaction diagram of a reinforced concrete column cross-section
% corresponding to a given axial force. The bisection root method is used.
% 
% OUTPUT: root:                 is a vector containing the neutral axis 
%                               depth, axial resistant force and bending
%                               resistance of a reinforced column
%                               cross-section as [c,FR,MR]
%
% INPUT:  cUno,cDos:            are the initial values of the neutral axis
%                               to commence iterations
%
%         fr:                   is the axial force resistance corresponding
%                               to the bending moment resistance for which
%                               the equilibrium condition sum F=0 is 
%                               established to extract its corresponding 
%                               bending moment resistance and neutral axis
%                               depth from the interaction diagram
%
%         E:                    Elasticity modulus of steel (4200 Kg/cm^2)
%
%         diam:                 is the cross-section diameter
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta1:                is determined as stablished by code (see
%                               Documentation)
%
%         ea:                   is the approximation error to terminate the
%                               root bisection method
%
%         av:                   ist the cross-section area of each rebar
%
%         rebarDisposition:     are the local coordinates of rebars over 
%                               the cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%

[eMecVar]=eleMecanicosRebarColsCirc(cUno,fdpc,diam,rec,E,rebarDisposition,...
                                    av,beta1);
                              
frt=eMecVar(1,1)+eMecVar(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eMecVar]=eleMecanicosRebarColsCirc(cDos,fdpc,diam,rec,E,rebarDisposition,...
                                    av,beta1);
frt=eMecVar(1,1)+eMecVar(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
            
[eMecVar]=eleMecanicosRebarColsCirc(c,fdpc,diam,rec,E,rebarDisposition,...
                                    av,beta1);
frt=eMecVar(1,1)+eMecVar(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%% beginning of loop %%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)

    if((raizUno*raizc)<0)
        cDos=c;
        [eMecVar]=eleMecanicosRebarColsCirc(cDos,fdpc,diam,rec,E,...
            rebarDisposition,av,beta1);
        frt=eMecVar(1,1)+eMecVar(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        [eMecVar]=eleMecanicosRebarColsCirc(cUno,fdpc,diam,rec,E,...
            rebarDisposition,av,beta1);
        
        frt=eMecVar(1,1)+eMecVar(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;

    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c==0
        c=0.001;
    end
    [eMecVar]=eleMecanicosRebarColsCirc(c,fdpc,diam,rec,E,rebarDisposition,...
                                    av,beta1);
    frt=eMecVar(1,1)+eMecVar(2,1);
    raizc=fr-frt;

    if (c~=0)
        es=abs((c-cu)/c);
    end

    if itdos>100 || ituno>100
        break;
    end
end
mrt=eMecVar(1,2)+eMecVar(2,2);
root=[c,frt,mrt];   
