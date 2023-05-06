function [root]=bisectionMrVarAsymRot(cUno,cDos,fr,E,h,b,fdpc,beta1,...
        ea,nv,number_rebars_sup,number_rebars_inf,number_rebars_izq,...
        number_rebars_der,rebarAvailable,op1,op2,op3,op4,...
        dispositionRebar,cp,RotCornerSec,gamma)

%------------------------------------------------------------------------
% Syntax:
% [raiz]=bisectionMrVarAsymRot(cUno,cDos,fr,E,h,b,fdpc,beta,...
%       ea,nv,number_rebars_sup,number_rebars_inf,number_rebars_izq,...
%       number_rebars_der,varDisponibles,op1,op2,op3,op4,...
%       disposicion_varillado,cp,RotCornerSec,gamma)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the neutral axis depth, axial and bending resistance
% from the interaction diagram of a rotated reinforced concrete column 
% cross-section with asymmetric reinforcement, for each for its points 
% with the aid of the bisection root method.
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
%                               depth from the interaction diagram. For
%                               columns, the value of this variable may go
%                               from the max resistance in compression of
%                               the cross-section (poc) to the max
%                               resistance in tension (pot)
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
%                               the rotated cross-section
%
%         rebarAvailable:       data base of commercial available rebars. 
%                               An array of size n# x 2 by default; in 
%                               format:
%                               [eight-of-an-inch, diam]
%
%         op1,op2,op3,op4:      eight-of-an-inch rebar to be placed for 
%                               each of the four boundaries of the 
%                               cross-section (upper boundary, lower 
%                               boundary, left boundary and right boundary)
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
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%% f(l) %%%%%%%%%%%%%%%%%%%%%

[eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
 number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,....
 op1,op2,op3,op4,b,h,cUno,fdpc,E,beta1,cp,RotCornerSec,gamma);
                              
frt=eMecVar(1,1)+eMecVar(2,1);
raizUno=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(u) %%%%%%%%%%%%%%%%%%%%%%% 
[eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
 number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,...
 op1,op2,op3,op4,b,h,cDos,fdpc,E,beta1,cp,RotCornerSec,gamma);

frt=eMecVar(1,1)+eMecVar(2,1);
raizDos=fr-frt;

%%%%%%%%%%%%%%%%%%%%%% f(xr) %%%%%%%%%%%%%%%%%%%%%%
c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));

if c==0
    c=1e-6;
end

[eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
 number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,...
 op1,op2,op3,op4,b,h,c,fdpc,E,beta1,cp,RotCornerSec,gamma);

frt=eMecVar(1,1)+eMecVar(2,1);
raizc=fr-frt;

%% Main loop
ituno=0;
itdos=0;

cu=cDos;
es=abs((c-cu)/c);
while(es>ea)
    
    if((raizUno*raizc)<0)
        cDos=c;
        [eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
         number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,...
         op1,op2,op3,op4,b,h,cDos,fdpc,E,beta1,cp,RotCornerSec,gamma);

        frt=eMecVar(1,1)+eMecVar(2,1);
        raizDos=fr-frt;

        ituno=ituno+1;

    elseif((raizUno*raizc)>0)
        cUno=c;
        [eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
         number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,...
         op1,op2,op3,op4,b,h,cUno,fdpc,E,beta1,cp,RotCornerSec,gamma);

        frt=eMecVar(1,1)+eMecVar(2,1);
        raizUno=fr-frt;

        itdos=itdos+1;
    end

    cu=c;

    c=cDos-(raizDos*(cUno-cDos)/(raizUno-raizDos));
    if c==0
        c=0.00001;
    end
    [eMecVar]=eleMecBarAsymRecRot(dispositionRebar,nv,number_rebars_sup,...
     number_rebars_inf,number_rebars_izq,number_rebars_der,rebarAvailable,...
     op1,op2,op3,op4,b,h,c,fdpc,E,beta1,cp,RotCornerSec,gamma);

    frt=eMecVar(1,1)+eMecVar(2,1);
    raizc=fr-frt;

    es=abs((c-cu)/c);

    if itdos>100 || ituno>100
        break;
    end
end
mrt=eMecVar(1,2)+eMecVar(2,2);
root=[c,frt,mrt];   
