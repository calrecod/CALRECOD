function [Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(t,...
        dimensionesColumna,rec,fy,npdiag,conditions,fdpc,E,beta1) 

%------------------------------------------------------------------------
% Syntax:
% [Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(t,...
%        dimensionesColumna,rec,fy,npdiag,conditions,fdpc,E,beta) 
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
% 
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram of an ISR reinforced column
% cross-section and its structural efficiency given some load conditions.
% 
% OUTPUT: Eft:                 Structural efficiency 
%
%         diagramaInteraccion: Interaction diagram coordinates for both 
%                              cross-section axis directions
%
%         tablaEficiencias:    is the resume table of results consisting of
%                              nloads rows and eight columns as:
%                              [Pu,Mux,Muy,PRx,PRy,MRx,MRy,Eff] 
%
%         cxy:                 is a vector containing the neutral axis depth
%                              of each cross-section direction according to
%                              the most critical load condition as [cx,cy]
%
% INPUT:  dimensionesColumna: cross-section dimensions of column (width 
%                             and height)
%
%         npdiag:             number of points to be analysed from the 
%                             interaction diagram
%
%         t:                  is the ISR width of a 1t-ISR for column 
%                             cross-sections
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
%         E:                  Elasticity modulus in units (Kg/cm^2)
%
%         c:                  neutral axis value (cm)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

ea=0.001; % approximation error for each point of the interaction diagram
if npdiag<3
    disp('Error: the number of points for the Interaction Diagram')
    disp('must be at least 3');
    return;
end
c_vector=zeros(npdiag,2);
diagramaInteraccion=zeros(npdiag,9);

h=dimensionesColumna(2);
b=dimensionesColumna(1);

rec_columna=rec;
act=(2*(b-2*rec(1))+2*(h-2*rec(2)))*t;

poc=(act*fy+fdpc*(b*h-act));
pot=act*fy;
df=(poc+pot)/(npdiag-1);
dfi=0.1/(npdiag-1);
c_vector(1,:)=4*h;

diagramaInteraccion(1,1)=-poc;
diagramaInteraccion(npdiag,1)=pot;
diagramaInteraccion(1,2)=0;
diagramaInteraccion(1,6)=0;
for sentido=1:2
    if (sentido==2)
        h=dimensionesColumna(1);
        b=dimensionesColumna(2);

        rec(1)=rec_columna(2);
        rec(2)=rec_columna(1);
    end
    diagramaInteraccion(1,4*sentido-1)=0.65*-poc;
    diagramaInteraccion(1,4*sentido)=0;
    for i=1:npdiag-1
        diagramaInteraccion(i+1,1)=-poc+i*df;
        cUno=0.001; 
        cDos=4*h;
        fr=diagramaInteraccion(i+1,1);
        
        [raiz]=bisectionMr4t(cUno,cDos,fr,E,t,t,t,t,h,b,rec,fdpc,beta1,ea);
        
        % Another option for computing the interaction diagram curves would
        % be using an analytical mathematical approach with the function:
        % 
        %[raiz]=bisectionMrAnalytISR(cUno,cDos,fr,E,t,h,b,rec,fdpc,beta1,ea);
        % 
        diagramaInteraccion(i+1,4*sentido-2)=raiz(3);
        c_vector(i+1,sentido)=raiz(1);
        
        %%%%%%%%%%%%%%% Reduced diagrams %%%%%%%%%%%%%%%%%%%%
        diagramaInteraccion(i+1,4*sentido-1)=(0.65+i*dfi)*diagramaInteraccion(i+1,1);
        diagramaInteraccion(i+1,4*sentido)=(0.65+i*dfi)*diagramaInteraccion(i+1,4*sentido-2);
                                            
        %%%%%%%%%%%%%%%%%%%%%%%%%% Eccentricities %%%%%%%%%%%%%%%%%%%%%%%%%%%

        diagramaInteraccion(i+1,4*sentido+1)=diagramaInteraccion(i+1,4*sentido)/...
                                         diagramaInteraccion(i+1,4*sentido-1);
    end

    diagramaInteraccion(npdiag,4*sentido-1)=0.75*pot;
    
    diagramaInteraccion(npdiag,4*sentido)=0;
    diagramaInteraccion(npdiag,2)=0;
    diagramaInteraccion(npdiag,6)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Efficiencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions(:,3:4)=abs(conditions(:,3:4)); % only with positive moments
                                          % is enough for symmetrical
                                          % designs, given that the ISR
                                          % considered in this function
                                          % is symmetrical ------------
[maxef,tablaEficiencias,cxy]=effRecColsLinearSearch(diagramaInteraccion,...
                                    conditions,pot,poc,c_vector);
                            
Eft=maxef;
