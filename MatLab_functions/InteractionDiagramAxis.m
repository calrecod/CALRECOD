function [diagramIntAxis1,pot,poc,cvectorX,newdispositionRebar,...
    newCoordCorners,newCP,gamma]=InteractionDiagramAxis...
    (npdiag,comborebar,b,h,fy,fdpc,beta1,E,number_rebars_sup,...
    number_rebars_inf,number_rebars_left,number_rebars_right,...
    rebarAvailable,dispositionRebar,Mux,Muy)

%------------------------------------------------------------------------
% Syntax:
% [diagramIntAxis1,pot,poc,cvectorX,newdispositionRebar,...
%  newCoordCorners,newCP,gamma]=InteractionDiagramAxis...
%  (npdiag,comborebar,b,h,fy,fdpc,beta1,E,number_rebars_sup,...
%  number_rebars_inf,number_rebars_left,number_rebars_right,...
%  rebarAvailable,dispositionRebar,Mux,Muy)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram with respect to a rotated
%          axis of a rectangular cross-section asymmetrically reinforced. 
% 
% OUTPUT: diagramIntAxis1:      is the array containing the interaction 
%                               diagram data whose rebars in tension are at
%                               the bottom of the cross-section. See
%                               Documentation
%
%         pot,poc:              are the max axial load resistance in
%                               tension and compression, respectively
%
%         newCP:                are the Plastic Center depth over the
%                               cross-section (with respect to the upper 
%                               outer most concrete fibre) in both axis
%                               directions. See documentation
% 
%         cvectorX:             is the neutral axis depth values along the
%                               interaction diagram of the axis in quest
%                               (from the upper cross-section corner
%                               downwards)
%
%         newdispositionRebar:  is the array containing the local
%                               coordinates of the rebars distributed over
%                               the rotated rectangular cross-section
%
%         gamma:                is the angle of rotation for the
%                               cross-section
%
% INPUT:  b,h:                  cross-section initial dimensions
%
%         rebarAvailable:       rebar database consisting of an array of 
%                               size [7,2] by default in format: 
%                               [#rebar,diam]
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
%         fdpc:                 is the reduced concrete's compressive
%                               strength by the application of the
%                               reduction resistance factor of 0.85
%                               (according to the ACI 318 code)
%
%         npdiag:               number of points to be computed for the
%                               definition of the interaction diagram
%
%         dispositionRebar:     is the vector containing the rebar local 
%                               coordinates over the cross-section [x,y] 
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-31
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
if npdiag<3
    disp('Error: the number of points for any axis of the Interaction') 
    disp('Surface must be at least 3');
    return;
end

op1=comborebar(1);
op2=comborebar(2);
op3=comborebar(3);
op4=comborebar(4);

diagramIntAxis1=zeros(npdiag,5);
cvectorX=zeros(npdiag,2);

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

nv=nv1+nv2+nv3+nv4;

arregloVar(1:nv1)=op1;
arregloVar(nv1+1:nv1+nv2)=op2;
arregloVar(nv1+nv2+1:nv1+nv2+nv3)=op3;
arregloVar(nv1+nv2+nv3+1:nv1+nv2+nv3+nv4)=op4;

%% Plastic Center in original coordinates of the reinforced cross-section
rebar=[dispositionRebar(:,1) dispositionRebar(:,2)];
for sentido=1:2
    if (sentido==2)
       h=dimensionesColumna(1);
       b=dimensionesColumna(2);
       dispositionRebar(:,1)=-rebar(:,2);
       dispositionRebar(:,2)=rebar(:,1);
    end
    % Plastic Center 
    [cp]=PlastiCenterAxis(fy,fdpc,b,h,dispositionRebar,arregloVar,...
                            rebarAvailable);
    CPaxis(1,sentido)=cp; % [cp_y,cp_x]
end

%% Rotation of the reinforced cross-section
h=dimensionesColumna(2);
b=dimensionesColumna(1);

alpha=rad2deg(atan(Mux/Muy));
if Muy<=0 && Mux>=0 || Muy<=0 && Mux<=0
    gamma=(90-alpha)+180;
elseif Muy>=0 && Mux<=0 || Muy>=0 && Mux>=0
    gamma=90-alpha; % This is the angle for the section to be rotated at
                    % so that the resultant moment of Mx and My is now Mx'
end

[newdispositionRebar,newCoordCorners,newCP]=rotReCol2(Mux,Muy,...
                                            rebar,b,h,CPaxis);
                                
%% Interaction diagram - Direction X'
    
poc=(act*fy+fdpc*(b*h-act));
pot=act*fy;
df=(poc+pot)/(npdiag-1);
dfi=0.1/(npdiag-1);

diagramIntAxis1(1,1)=-poc;
diagramIntAxis1(1,3)=0.65*-poc;
diagramIntAxis1(1,2)=1e-7;
diagramIntAxis1(1,4)=1e-7;

cvectorX(1,:)=4*max(newCoordCorners(:,2));
cvectorX(npdiag,:)=0;

ea=0.001;
for i=1:npdiag-1
    diagramIntAxis1(i+1,1)=-poc+i*df;
    cUno=0.001; 
    cDos=4*max(newCoordCorners(:,2));
    fr=diagramIntAxis1(i+1,1);

    [raiz]=bisectionMrVarAsymRot(cUno,cDos,fr,E,h,b,fdpc,beta1,ea,nv,...
        number_rebars_sup,number_rebars_inf,number_rebars_left,...
        number_rebars_right,rebarAvailable,op1,op2,op3,op4,...
        newdispositionRebar,newCP(1),newCoordCorners,gamma);

    diagramIntAxis1(i+1,2)=raiz(3);
    cvectorX(i+1,1)=raiz(1);

    %%%%%%%%%%%%%%% Reduced diagramas %%%%%%%%%%%%%%%%%%%%
    diagramIntAxis1(i+1,3)=(0.65+i*dfi)*diagramIntAxis1(i+1,1);
    diagramIntAxis1(i+1,4)=(0.65+i*dfi)*diagramIntAxis1(i+1,2);

    %%%%%%%%%%%%%%%%%%%%%%%%% Eccentricities %%%%%%%%%%%%%%%%%%%%%%%%%%

    diagramIntAxis1(i+1,5)=diagramIntAxis1(i+1,4)/...
                           diagramIntAxis1(i+1,3);
end
diagramIntAxis1(npdiag,1)=pot;
diagramIntAxis1(npdiag,3)=0.75*pot;
diagramIntAxis1(npdiag,4)=1e-7;
diagramIntAxis1(npdiag,2)=1e-7;
