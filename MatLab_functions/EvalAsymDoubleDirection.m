function [diagramaInteraccion1,diagramaInteraccion2,pot,poc,cp_axis,cvector1,...
    cvector2]=EvalAsymDoubleDirection(npdiag,comborebar,b,h,...
    fy,fdpc,beta,E,number_rebars_sup,number_rebars_inf,number_rebars_left,...
    number_rebars_right,rebarAvailable,dispositionRebar,concreteCover)

%------------------------------------------------------------------------
% Syntax:
% [diagramaInteraccion1,diagramaInteraccion2,poc,pot,cp_axis,cvector1,...
% cvector2]=EvalAsymDoubleDirection(npdiag,comborebar,b,h,fy,fdpc,beta,...
% E,number_rebars_sup,number_rebars_inf,number_rebars_left,
% number_rebars_right,rebarAvailable,dispositionRebar,concreteCover)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagrams of an asymmetrically 
% reinforced cross-section design option, in both directions (positive and
% negative) of benidng moments.
% 
% OUTPUT: diagramaInteraccion1:  is the array containing that interaction 
%                               diagram data for positive bending moments
%
%         diagramaInteraccion2:  is the array containing that interaction 
%                               diagram data for negative bending moments
%
%         pot,poc,:             are the max axial load resistances for the
%                               reinforced cross-section in tension and 
%                               compression, respectively
%
%         cp_axis:              are the Plastic Center depth over the
%                               cross-section (with respect to the upper 
%                               outer most concrete fibre) in both axis
%                               directions. See documentation
% 
%         cvector1, cvector2:   are the vectors containing the neutral axis
%                               depth values along each interaction diagram 
%                               of both cross-section's axis
%
% INPUT:  b,h:                  cross-section initial dimensions
%
%         rebarAvailable:       rebar database consisting of an array of 
%                               size [7,3] by default in format: 
%                               [#rebar,diam,unit-weight]
%
%         number_rebars_sup,
%         number_rebars_inf,
%         number_rebars_left,
%         number_rebars_right:  number of rebars to placed on each of the 
%                               cross-section boundaries
%
%         comborebar:           vector of size [1,4] which are the
%                               types of rebar for each of the 
%                               four cross-section boundaries (upper 
%                               boundary, lower boundary, left side and 
%                               right side, respectively)
%
%         npdiag:               number of points to be computed for the
%                               definition of the interaction diagram
%
%         concreteCover:        concrete cover along each cross-section's
%                               axis (x,y)
%
%         dispositionRebar:     is the vector containing the rebar local 
%                               coordinates over the cross-section [x,y] 
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
if npdiag<3
    disp('Error: the number of points for the Interaction Diagram must be 3 or higher');
    return;
end
cp_axis=[];

op1=comborebar(1);
op2=comborebar(2);
op3=comborebar(3);
op4=comborebar(4);

diagramaInteraccion1=zeros(npdiag,9);
cvector1=zeros(npdiag,2);

dimensionesColumna=[b h];

h=dimensionesColumna(2);
b=dimensionesColumna(1);

act=number_rebars_sup*(rebarAvailable(op1,2)^2*pi/4)+...
    number_rebars_inf*(rebarAvailable(op2,2)^2*pi/4)+...
    number_rebars_left*(rebarAvailable(op3,2)^2*pi/4)+...
    number_rebars_right*(rebarAvailable(op4,2)^2*pi/4);

nv1=number_rebars_sup;
nv2=number_rebars_inf;
nv3=number_rebars_left;
nv4=number_rebars_right;

nb4=[nv1,nv2,nv4,nv4];

nv=nv1+nv2+nv3+nv4;

arregloVar=zeros(nv,1);
for i=1:nv1
    arregloVar(i)=op1;
end

for i=1:nv2
    arregloVar(nv1+i)=op2;
end

for i=1:nv3
    arregloVar(nv1+nv2+i)=op3;
end

for i=1:nv4
    arregloVar(nv1+nv2+nv3+i)=op4;
end

rebar=[dispositionRebar(:,1) dispositionRebar(:,2)];

%%%%%%%%%%%%%%%% INTERACTION DIAGRAMS-DIRECTION 1 %%%%%%%%%%%%%%%%%%%%%
    
poc=(act*fy+fdpc*(b*h-act));
pot=act*fy;
df=(poc+pot)/(npdiag-1);
dfi=0.1/(npdiag-1);

diagramaInteraccion1(1,1)=-poc;
diagramaInteraccion1(npdiag,1)=pot;
diagramaInteraccion1(1,2)=1e-7;
diagramaInteraccion1(1,6)=1e-7;

cvector1(1,:)=4*h;
cvector1(npdiag,:)=0;

ea=0.001;

for sentido=1:2
    if (sentido==2)
       h=dimensionesColumna(1);
       b=dimensionesColumna(2);
       dispositionRebar(:,1)=rebar(:,2); % the concrete section is rotated
       dispositionRebar(:,2)=-rebar(:,1);% counter clock-wise, but not 
                                         % the rebars
    end
    % Plastic Center 
    [cp]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,arregloVar,...
                            rebarAvailable);
    cp_axis=[cp_axis,cp];
   
    diagramaInteraccion1(1,4*sentido-1)=0.65*-poc;
    diagramaInteraccion1(1,4*sentido)=1e-7;
    for i=1:npdiag-1
        diagramaInteraccion1(i+1,1)=-poc+i*df;
        cUno=0.001; 
        cDos=4*h;
        fr=diagramaInteraccion1(i+1,1);

        [raiz]=bisectionMrVarAsymm(cUno,cDos,fr,E,h,b,fdpc,beta,ea,nv,...
            number_rebars_sup,number_rebars_inf,number_rebars_left,...
            number_rebars_right,rebarAvailable,op1,op2,op3,op4,...
            dispositionRebar,cp);

        diagramaInteraccion1(i+1,4*sentido-2)=raiz(3);
        cvector1(i+1,sentido)=raiz(1);

        %%%%%%%%%%%%%%% Reduced diagramas %%%%%%%%%%%%%%%%%%%%
        diagramaInteraccion1(i+1,4*sentido-1)=(0.65+i*dfi)*diagramaInteraccion1(i+1,1);
        diagramaInteraccion1(i+1,4*sentido)=(0.65+i*dfi)*diagramaInteraccion1(i+1,4*sentido-2);
                                                          
        %%%%%%%%%%%%%%% Eccentricities %%%%%%%%%%%%%%%%%%%%%%%%%%

        diagramaInteraccion1(i+1,4*sentido+1)=diagramaInteraccion1(i+1,4*sentido)/...
                                         diagramaInteraccion1(i+1,4*sentido-1);
    end

    diagramaInteraccion1(npdiag,4*sentido-1)=0.75*pot;
    diagramaInteraccion1(npdiag,4*sentido)=1e-7;
    diagramaInteraccion1(npdiag,4*sentido-2)=1e-7;
    diagramaInteraccion1(npdiag,4*sentido-2)=1e-7;
end

%%%%%%%%%%%%%%%% INTERACTION DIAGRAMS - DIRECTION 2 %%%%%%%%%%%%%%%%%%%

numberRebars1=number_rebars_inf;
numberRebars2=number_rebars_sup;
numberRebars3=number_rebars_right;
numberRebars4=number_rebars_left;

RebarTypeIndex1=op2;
RebarTypeIndex2=op1;
RebarTypeIndex3=op4;
RebarTypeIndex4=op3;
[dispositionRebar,separacion_hor1,separacion_hor2,...
            separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
            h,concreteCover,nv,numberRebars1,numberRebars2,...
            numberRebars3,numberRebars4,rebarAvailable,RebarTypeIndex1,...
            RebarTypeIndex2,RebarTypeIndex3,RebarTypeIndex4);
        
arregloVar=zeros(nv,1);

nv1=numberRebars1;
nv2=numberRebars2;
nv3=numberRebars3;
nv4=numberRebars4;
for i=1:nv1
    arregloVar(i)=RebarTypeIndex1;
end

for i=1:nv2
    arregloVar(nv1+i)=RebarTypeIndex2;
end

for i=1:nv3
    arregloVar(nv1+nv2+i)=RebarTypeIndex3;
end

for i=1:nv4
    arregloVar(nv1+nv2+nv3+i)=RebarTypeIndex4;
end

diagramaInteraccion2=zeros(npdiag,9);
cvector2=zeros(npdiag,2);

diagramaInteraccion2(1,1)=-poc;
diagramaInteraccion2(npdiag,1)=pot;
diagramaInteraccion2(1,2)=-1e-7;
diagramaInteraccion2(1,6)=-1e-7;

cvector2(1,:)=4*h;
cvector2(npdiag,:)=0;

h=dimensionesColumna(2);
b=dimensionesColumna(1);

rebar=[dispositionRebar(:,1),dispositionRebar(:,2)];

for sentido=1:2
    if (sentido==2)
       h=dimensionesColumna(1);
       b=dimensionesColumna(2);
       dispositionRebar(:,1)=rebar(:,2);
       dispositionRebar(:,2)=-rebar(:,1);
    end
    % Plastic Center
    [cp]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,arregloVar,...
                            rebarAvailable);
                        
    diagramaInteraccion2(1,4*sentido-1)=0.65*-poc;
    diagramaInteraccion2(1,4*sentido)=-1e-7;
    for i=1:npdiag-1
        diagramaInteraccion2(i+1,1)=-poc+i*df;
        cUno=0.001; 
        cDos=4*h;
        fr=diagramaInteraccion2(i+1,1);
        
        [raiz]=bisectionMrVarAsymm(cUno,cDos,fr,E,h,b,fdpc,beta,ea,nv,...
            numberRebars1,numberRebars2,numberRebars3,...
            numberRebars4,rebarAvailable,RebarTypeIndex1,RebarTypeIndex2,...
            RebarTypeIndex3,RebarTypeIndex4,dispositionRebar,cp);

        diagramaInteraccion2(i+1,4*sentido-2)=-raiz(3);
        cvector2(i+1,sentido)=raiz(1);

        %%%%%%%%%%%%%%% Reduced diagramas %%%%%%%%%%%%%%%%%%%%
        diagramaInteraccion2(i+1,4*sentido-1)=(0.65+i*dfi)*diagramaInteraccion2(i+1,1);
        diagramaInteraccion2(i+1,4*sentido)=(0.65+i*dfi)*diagramaInteraccion2(i+1,4*sentido-2);
                                                          
        %%%%%%%%%%%%%%%%% Eccentricities %%%%%%%%%%%%%%%%%%%%%%%%%%
        diagramaInteraccion2(i+1,4*sentido+1)=diagramaInteraccion2(i+1,4*sentido)/...
                                         diagramaInteraccion2(i+1,4*sentido-1);
    end

    diagramaInteraccion2(npdiag,4*sentido-1)=0.75*pot;
    diagramaInteraccion2(npdiag,4*sentido)=-1e-7;
    diagramaInteraccion2(npdiag,4*sentido-2)=-1e-7;
    diagramaInteraccion2(npdiag,4*sentido-2)=-1e-7;
end

