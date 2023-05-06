function [maxef,diagramaInteraccion,tablaEficiencias,c]=diagCircColsRebar(ast,...
        dispositionRebar,diam,rec,fy,npdiag,conditions,fdpc,av,E,beta1) 

%------------------------------------------------------------------------
% Syntax:
% [maxef,diagramaInteraccion,tablaEficiencias,c]=diagCircColsRebar(ast,...
%  dispositionRebar,diam,rec,fy,npdiag,conditions,fdpc,av,E,fdpc,beta) 
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram of a rebar design on a 
% circular column  and its structural efficiency given some 
% load conditions.
% 
% OUTPUT: maxef:               Structural efficiency 
%
%         diagramaInteraccion: Interaction diagram coordinates 
%
%         tablaEficiencias:    is the resume table of results consisting of
%                              nloads rows and five columns as:
%                              [Pu,Mu,PR,MR,Eff] 
%
%         c:                   is the neutral axis depth corresponding 
%                              to the most critical load condition
%
% INPUT:  diam:               cross-section diameter 
%
%         npdiag:             number of points to be analysed from the 
%                             interaction diagram
%
%         rec:                is the concrete cover
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code
%
%         beta1:              is determined as established in ACI 318 code
%                             (see Documentation)
%
%         E:                  Modulus of Elasticity of the reinforcement
%                             steel
%
%         conditions:         is the array containing the load conditions
%                             - size: nloads x 3, in format [nload,Pu,Mu]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

if npdiag<3
    disp('Error: the number of points for the Interaction Diagram must be 3 or higher');
    return;
end
c_vector=zeros(npdiag,1);
diagramaInteraccion=zeros(npdiag,5);

acol=0.25*pi*diam^2;

poc=(ast*fy+fdpc*(acol-ast));
pot=ast*fy;
df=(poc+pot)/(npdiag-1);
dfi=0.1/(npdiag-1);
c_vector(1,1)=4*diam;

diagramaInteraccion(1,1)=-poc;
diagramaInteraccion(npdiag,1)=pot;
diagramaInteraccion(1,2)=0;

diagramaInteraccion(1,3)=0.65*-poc;
diagramaInteraccion(1,4)=0;
for i=1:npdiag-1
    diagramaInteraccion(i+1,1)=-poc+i*df;
    cUno=0.001; 
    cDos=diam;
    fr=diagramaInteraccion(i+1,1);
    [raiz]=bisectionMrRebarColsCir(cUno,cDos,fr,E,diam,fdpc,beta1,...
                             0.001,av,dispositionRebar,rec);  

    diagramaInteraccion(i+1,2)=raiz(3);
    c_vector(i+1,1)=raiz(1);

    %%%%%%%%%%%%%%%%%%%%%%%% Reduced diagrams %%%%%%%%%%%%%%%%%%%%%%%%
    diagramaInteraccion(i+1,3)=(0.65+i*dfi)*diagramaInteraccion(i+1,1);
    diagramaInteraccion(i+1,4)=(0.65+i*dfi)*diagramaInteraccion(i+1,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%% Eccentricity %%%%%%%%%%%%%%%%%%%%%%%%%%

    diagramaInteraccion(i+1,5)=diagramaInteraccion(i+1,4)/...
                               diagramaInteraccion(i+1,3);
end

diagramaInteraccion(npdiag,3)=0.75*pot;

diagramaInteraccion(npdiag,4)=0;
diagramaInteraccion(npdiag,2)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Eficiencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxef,tablaEficiencias,c]=effCircColsLS(diagramaInteraccion,...
                                        conditions,c_vector);
