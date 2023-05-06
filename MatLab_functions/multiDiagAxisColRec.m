function [tablaEff,iloadmax,maxLoadCondition,maxgamma,diagramIntAxis1,...
    rotdispositionRebar,rotSection,cmax,CP]=multiDiagAxisColRec(b,h,...
    load_conditions,comborebar,npdiag,fy,fdpc,beta1,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
    dispositionRebar)

%------------------------------------------------------------------------
% Syntax:
% [tablaEff,iloadmax,maxLoadCondition,maxgamma,diagramIntAxis1,...
%  rotdispositionRebar,rotSection,cmax,CP]=multiDiagAxisColRec(b,h,...
%  load_conditions,comborebar,npuntos,fy,fdpc,beta1,E,numberRebars1,...
%  numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
%  dispositionRebar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To determine the structural efficiency of the most critical of
% the given biaxial load conditions for a rectangular asymmetrically
% reinforced cross-section column. A 3D analysis is carried on for the
% computation of the interaction diagrams in the plane in which each of the
% load conditions are applied.
% 
% INPUT:  b,h:                  are the cross-section dimensions: width and
%                               heigth
%
%         load_conditions:      load conditions vector of size nloads x 4
%                               in format: [nload,Pu,Mux,Muy]
%
%         comborebar:           is the vector containing the four rebar
%                               diamemters' indices (from the given rebar 
%                               database table), in format:
%
%   [index-bar-upper, index-bar-lower, index-bar-left,index-bar-right]
%
%         npdiag:               number of points to be computed for the
%                               interaction diagram 
%
%         fy:                   list of types of rebar distributed over the
%                               cross-section
%
%         fdpc:                 is the concrete compressive strength 
%
%         beta1:                is determined as prescribed in the 
%                               ACI 318-19 code (see Documentation)
%                               
%
%         E:                    is the Modulus of Elasticity of the
%                               reinforcing steel
%
%         numberRebars1,
%         numberRebars2,
%         numberRebars3,
%         numberRebars4:        are the number of rebars on each
%                               cross-section's boundary (upper boundary,
%                               lower boundary, left boundary, right
%                               boundary)
%
%         rebarAvailable:       is the rebar database table
%
%         dispositionRebar:     is the array that contains the rebar
%                               coordinates over the cross-section
%
% OUTPUT: tablaEff:             is the resume table of the structural
%                               efficiency analysis for each load
%                               condition. Size: nloads x 5, in format:
%                               [Pu,Mu,Pr,Mr,Eff] 
%
%         iloadmax:             is the index of the most critical load
%                               condition (according to the order in which
%                               they were given in the "load condition"
%                               array
%
%         maxLoadCondition:     is the structural efficiency analysis of
%                               the most critical load condition. Size: 
%                               1 x 3, in format: [1,Pu,Mxy] where:
%                               Mxy=sqrt(Mux^2+Muy^2)
%
%         maxgamma:             is the angle of rotation that the
%                               cross-section had to undergo for the most 
%                               critical load condition 
%
%         diagramIntAxis1:      is the array containing the interaction
%                               diagram's data for the most critical of the
%                               given load conditions
%
%         rotdispositionRebar:  is the array containing the rotated rebar
%                               coordinates over the cross-section 
%                               corresponding to the most critical of the 
%                               load conditions
%
%         rotSection:           are the coordinates of the rotated
%                               cross-section for the most critical of the
%                               given load conditions
%
%         cmax:                 is the neutral axis depth value of the
%                               rotated cross-section for the most critical
%                               of the load conditions
%
%         CP:                   are the Plastic Center depth values
%                               corresponding to the rotated cross-section
%                               for the most critical of the given load
%                               conditions
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-04-12
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

maxef=-100000;
nloads=length(load_conditions(:,1));
for i=1:nloads
    %% Interaction diagram in the load direction
    Pu=load_conditions(i,2);
    Mux=load_conditions(i,3);
    Muy=load_conditions(i,4);

    [diagramIntAxis,pot,poc,cvector1,newdispositionRebar,...
    newCoordCorners,newCP,gamma]=InteractionDiagramAxis...
    (npdiag,comborebar,b,h,fy,fdpc,beta1,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,...
    rebarAvailable,dispositionRebar,Mux,Muy);

    newCoordCorners(5,:)=newCoordCorners(1,:);    % to close the drawing
                                                  % of the rotated section
    %% Structural efficiency
    Mxy=sqrt(Mux^2+Muy^2);
    trans_load_condition=[1 Pu Mxy];
    [ef,tablaEffLoad,c]=effColsRot1DirecLS(diagramIntAxis,...
                                trans_load_condition,cvector1);
    tablaEff(i,:)=tablaEffLoad;
    if ef>maxef
        maxgamma=gamma;
        maxef=ef;
        iloadmax=i;
        diagramIntAxis1=diagramIntAxis;
        rotdispositionRebar=newdispositionRebar;
        rotSection=newCoordCorners;
        cmax=c;
        CP=newCP;
        maxLoadCondition=trans_load_condition;
    end
end