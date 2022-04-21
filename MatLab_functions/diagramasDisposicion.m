function [diagrama,mexef,eficiencia,cxy]=diagramasDisposicion(As,b,h,E,npuntos,...
                       fdpc,nv,beta,ov,av,disposicion_varillado,load_conditions)
                            
%------------------------------------------------------------------------
% Syntax:
% [diagrama,mexef,eficiencia,cxy]=diagramasDisposicion(As,b,h,E,npuntos,...
%        fdpc,nv,beta,ov,av,disposicion_varillado,load_conditions)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram of a symmetric rebar option, 
% as well as the structural efficiency given certain load conditions.
% 
% OUTPUT: diagrama:             is the interaction diagram data
%
%         maxef:                is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         eficiencia:           is a table containing the structural 
%                               efficiency analysis data: size = [nload,8],
%                               in format: 
%           _____________________________________________________
%           [nload,Pu,Mux,Muy,P{Rx},M{Rx},P{Ry},M{Ry},efficiency]
%           _____________________________________________________
%
%         cxy:                  are the neutral axis depth values 
%                               corresponding to the critical load condition, 
%                               for both axis directions: [cx,cy]
%
% INPUT:  b,h:                  given cross-section dimensions
%
%         av,ov,nv:             are the rebar area and number of rebar of 
%                               the current rebar option, and the number of
%                               rebars, respectively
%
%         load_conditions:      is the array containing the load conditions:
%                               size = [nload,4] in format: [nload,Pu,Mux,Muy]
%
%         beta:                 is determined as specified in code ACI,
%                               according to the f'c used
%
%         fdpc:                 equal to 0.85f'c according to code
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------


diagrama=zeros(npuntos+2,9);
c_vector_bar=zeros(npuntos+2,2);

fy=4200;
Ac=b*h;
pot=As*fy*0.001;

poc=(fdpc*(Ac-As)+fy*As)*0.001;
df=(pot+poc)/(npuntos);

diagrama(1,1)=-poc;
diagrama(npuntos+2,1)=pot;
diagrama(1,2)=0;
diagrama(npuntos+2,2)=0;

c_vector_bar(1,1)=4*h;
c_vector_bar(1,2)=4*h;

ncondiciones=length(load_conditions(:,1));

varillado=[disposicion_varillado(:,1) disposicion_varillado(:,2)];
dimensionesColumna=[b h];
for sentido=1:2
   if (sentido==2)
       b=dimensionesColumna(2);
       h=dimensionesColumna(1);
       disposicion_varillado(:,1)=-varillado(:,2);
       disposicion_varillado(:,2)=varillado(:,1);
   end
for i=1:npuntos
    diagrama(i+1,1)=(-poc+i*df);
    fr=diagrama(i+1,1);
    cUno=0.001;
    cDos=4*h;

    [raiz]=bisectionMrSymRebarCols(cUno,cDos,fr,E,h,b,fdpc,beta,...
                             0.001,nv,ov,av,disposicion_varillado);  
    
    diagrama(i+1,4*sentido-2)=raiz(3);

    c_vector_bar(i+1,sentido)=raiz(1);
    
    %%%%%%%%%%%%%%% diagramas reducidos por FR %%%%%%%%%%%%%%%%%%%%
    
    if diagrama(i+1,4*sentido-2)>diagrama(i,4*sentido-2)
       diagrama(i+1,4*sentido-1)=0.65*diagrama(i+1,1);
       diagrama(i+1,4*sentido)=0.65*diagrama(i+1,4*sentido-2);
    else
        
        diagrama(i+1,4*sentido-1)=0.75*diagrama(i+1,1);
        diagrama(i+1,4*sentido)=0.75*diagrama(i+1,4*sentido-2);
    end 
%%%%%%%%%%%%%%%%%%%%%%%% Excentriciddes %%%%%%%%%%%%%%%%%%%%%%%%%%%

   diagrama(i+1,4*sentido+1)=diagrama(i+1,4*sentido)/...
                                     diagrama(i+1,4*sentido-1);
end
  
  diagrama(1,4*sentido-1)=0.65*-poc;
  diagrama(npuntos+2,4*sentido-1)=0.75*pot;
  diagrama(1,4*sentido)=0;
  diagrama(npuntos+2,4*sentido)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Eficiencias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mexef,eficiencia,cxy]=eficienciaRec(diagrama,load_conditions,pot,poc,c_vector_bar);
