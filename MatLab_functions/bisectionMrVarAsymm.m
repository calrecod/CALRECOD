function [raiz]=bisectionMrVarAsymm(cUno,cDos,fr,E,h,b,fdpc,beta,...
        ea,nv,nRebarsTop,nRebarsBot,nRebarsL,...
        nRebarsR,rebarAvailable,op1,op2,op3,op4,...
        dispositionRebar,cp)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrVarAsymm(cUno,cDos,fr,E,h,b,fdpc,beta,...
%       ea,nv,nRebarsTop,nRebarsBot,nRebarsL,...
%       nRebarsR,varDisponibles,op1,op2,op3,op4,...
%       dispositionRebar,cp)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth, axial and bending resistance
% from the interaction diagram of a reinforced concrete column cross-section
% with asymmetric reinforcement, for each for its points with the aid of 
% the bisection root method.
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
%                               bending moment resistance and neutral axis
%                               depth from the interaction diagram
%
%         E:                    Elasticity modulus of steel
%
%         b,h:                  cross-section dimensions
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta:                 is determined as stablished by code 
%                               (see Documentation)
%
%         ea:                   is the approximation error to terminate the
%                               root bisection method
%
%         nv:                   is the number of rebars to be placed over
%                               the cross-section
%
%         number_rebars_sup,
%         number_rebars_inf,
%         number_rebars_left,
%         number_rebars_right:  are the number of rebars to be placed for 
%                               each of the cross-section boundaries
%
%         dispositionRebar:     are the local coordinates of rebars over 
%                               the cross-section
%
%         rebarAvailable:       data base of commercial available rebars. 
%                               An array of size n# x 2 (by default); 
%                               in format [#rebar,diam,unit-weight]
%
%         op1,op2,op3,op4:      types of rebar to be placed for each of the
%                               four boundaries of the cross-section 
%                               (upper boundary, lower boundary, left 
%                               boundary and right boundary)
%
%         cp:                   Plastic Center location for each axis
%                               direction of the column cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%

[eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
        nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                      cUno,fdpc,E,beta,cp);
                              
frt=eMecVar(1,1)+eMecVar(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
        nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                      cDos,fdpc,E,beta,cp);
frt=eMecVar(1,1)+eMecVar(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));

if c==0
    c=0.00000001;
end
a=beta*c;

[eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
        nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                      c,fdpc,E,beta,cp);
frt=eMecVar(1,1)+eMecVar(2,1);
raizc=fr-frt;

%%%%%%%%%%%%%%%%%%%%%%%%%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)
    
    if((raizUno*raizc)<0)
        cDos=c;
        [eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
            nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                          cDos,fdpc,E,beta,cp);
        frt=eMecVar(1,1)+eMecVar(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        [eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
                nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                              cUno,fdpc,E,beta,cp);
            frt=eMecVar(1,1)+eMecVar(2,1);
            raizUno=fr-frt;

            itdos=itdos+1;
    end

    cu=c;

    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c==0
        c=0.00001;
    end
    [eMecVar]=eleMecanicosVarAsymm(dispositionRebar,nv,nRebarsTop,nRebarsBot,...
            nRebarsL,nRebarsR,rebarAvailable,op1,op2,op3,op4,b,h,...
                                          c,fdpc,E,beta,cp);
    frt=eMecVar(1,1)+eMecVar(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);

    if itdos>100 || ituno>100
        break;
    end
end
mrt=eMecVar(1,2)+eMecVar(2,2);
raiz=[c,frt,mrt];   
