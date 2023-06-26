function [diagramIntAxis1,pot,poc,cvectorX,newcoordISRdisc,...
    newCoordCorners,newCP,gamma]=InteracDiagAxisISR(npdiag,t,beta1,...
    dimensionsColumn,rec,fy,fdpc,E,Mux,Muy) 
%------------------------------------------------------------------------
% Syntax:
% [diagramIntAxis1,pot,poc,cvectorX,newcoordISRdisc,newCoordCorners,...
% newCP,gamma]=InteracDiagAxisISR(npdiag,t,beta1,dimensionsColumn,rec,...
% fy,fdpc,E,Mux,Muy)  
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
% 
%------------------------------------------------------------------------
% PURPOSE: To compute the interaction diagram of an ISR reinforced column
% cross-section and its structural efficiency given some load conditions.
% 
% Note: The cross-section is rotated according to the given pair of bending
% moments so that the interaction diagram corresponds to the axis of action
% of such pair of moments's axis.
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
% INPUT:  npdiag:             number of points to be analysed from the 
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
%         E:                  Modulus of Elasticity of reinforcing steel
%
%         Mux,Muy:            Bending moment components with respect to
%                             each cross-section local axis
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-22
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

ea=0.001; % approximation error for each point of the interaction diagram
if npdiag<3
    disp('Error: the number of points for the Interaction Diagram')
    disp('must be at least 3');
    return;
end

%% Computing the non-rotated coordinates of the discrete ISR
b=dimensionsColumn(1);
h=dimensionsColumn(2);
ndt=100;
[da4,coordISRdisc]=discreteISRColRec(t,b,h,rec,ndt);

%% Rotation of the reinforced cross-section
h=dimensionsColumn(2);
b=dimensionsColumn(1);

alpha=rad2deg(atan(Mux/Muy));
if Muy<=0 && Mux>=0 || Muy<=0 && Mux<=0
    gamma=(90-alpha)+180;
elseif Muy>=0 && Mux<=0 || Muy>=0 && Mux>=0
    gamma=90-alpha; % This is the angle for the section to be rotated at
                    % so that the resultant moment of Mx and My is now Mx'
end

[newcoordISRdisc,newCoordCorners,newCP]=rotReCol2(Mux,Muy,...
                           coordISRdisc,b,h,[0.5*h,0.5*b]);
                                        
%% Interaction diagram - Direction X'
cvectorX=zeros(npdiag,1);
diagramIntAxis1=zeros(npdiag,5);

act=(2*(b-2*rec(1))+2*(h-2*rec(2)))*t;

poc=(act*fy+fdpc*(b*h-act));
pot=act*fy;
df=(poc+pot)/(npdiag-1);
dfi=0.1/(npdiag-1);

cvectorX(1,:)=4*max(newCoordCorners(:,2));
cvectorX(npdiag,:)=0;

diagramIntAxis1(1,1)=-poc;
diagramIntAxis1(1,3)=0.65*-poc;
diagramIntAxis1(1,2)=1e-7;
diagramIntAxis1(1,4)=1e-7;
for i=1:npdiag-1
    diagramIntAxis1(i+1,1)=-poc+i*df;
    cUno=0.001; 
    cDos=4*max(newCoordCorners(:,2));
    fr=diagramIntAxis1(i+1,1);

    [raiz]=bisectionMrSymISRRot(cUno,cDos,fr,E,h,b,fdpc,beta1,ea,4*ndt,...
        da4,newcoordISRdisc,newCP(1),newCoordCorners,gamma);

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

