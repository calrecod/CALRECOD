function [Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(t,...
        dimensionesColumna,rec,fy,npuntos,conditions,fdpc,E,beta) 

%------------------------------------------------------------------------
% Syntax:
% [Eft,diagramaInteraccion,tablaEficiencias,cxy]=widthEfficiencyCols(t,...
%        dimensionesColumna,rec,fy,npuntos,conditions,fdpc,E,beta) 
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
%         npuntos:            number of points to be analysed from the 
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

ndiagrama=npuntos+1;
c_vector=zeros(ndiagrama,1);

ncondiciones=length(conditions(:,1));
diagramaInteraccion=zeros(ndiagrama,9);

h=dimensionesColumna(2);
b=dimensionesColumna(1);

rec_columna=rec;
act=(2*(b-2*rec(1))+2*(h-2*rec(2)))*t;

poc=(act*fy+fdpc*(b*h-act))*0.001;
pot=act*fy*0.001;
df=(poc+pot)/npuntos;

c_vector(1)=4*h;
diagramaInteraccion(1,1)=-poc;
diagramaInteraccion(npuntos+1,1)=pot;
diagramaInteraccion(1,2)=0;
diagramaInteraccion(npuntos+1,2)=0;
    
for sentido=1:2
    if (sentido==2)
        h=dimensionesColumna(1);
        b=dimensionesColumna(2);

        rec(1)=rec_columna(2);
        rec(2)=rec_columna(1);
    end
    for i=1:ndiagrama-1
        diagramaInteraccion(i+1,1)=-poc+i*df;
        cUno=0.001; 
        cDos=4*h;
        fr=diagramaInteraccion(i+1,1);
        [raiz]=bisectionMr4t(cUno,cDos,fr,E,t,t,t,t,h,b,rec,fdpc,beta,0.001);  

        if (i==0)
            diagramaInteraccion(i+1,4*sentido-2)=0.001;
        elseif(i==ndiagrama-1)
            diagramaInteraccion(i+1,4*sentido-2)=0.001;
        else
            diagramaInteraccion(i+1,4*sentido-2)=raiz(3);
        end
        c_vector(i+1)=raiz(1);
        
        %%%%%%%%%%%%%%% diagramas reducidos por FR %%%%%%%%%%%%%%%%%%%%
        if diagramaInteraccion(i+1,4*sentido-2)>diagramaInteraccion(i,4*sentido-2)
        %if fr<0   
           diagramaInteraccion(i+1,4*sentido-1)=0.65*diagramaInteraccion(i+1,1);
           diagramaInteraccion(i+1,4*sentido)=0.65*diagramaInteraccion(i+1,4*sentido-2);
        else

            diagramaInteraccion(i+1,4*sentido-1)=0.75*diagramaInteraccion(i+1,1);
            diagramaInteraccion(i+1,4*sentido)=0.75*diagramaInteraccion(i+1,4*sentido-2);
        end                                           
        %%%%%%%%%%%%%%%%%%%%%%%%%% Excentriciddes %%%%%%%%%%%%%%%%%%%%%%%%%%%

        diagramaInteraccion(i+1,4*sentido+1)=diagramaInteraccion(i+1,4*sentido)/...
                                         diagramaInteraccion(i+1,4*sentido-1);
    end

    diagramaInteraccion(1,4*sentido-1)=0.65*-poc;
    diagramaInteraccion(ndiagrama,4*sentido-1)=0.75*pot;
    diagramaInteraccion(1,4*sentido)=0;
    diagramaInteraccion(npuntos+1,4*sentido)=0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Eficiencias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[maxef,diagramaInteraccion,tablaEficiencias,cxy]=eficienciaRecISRCols(ncondiciones,...
    diagramaInteraccion,conditions,pot,poc,c_vector);

Eft=maxef;