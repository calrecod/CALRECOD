function [maxef,tablaEff,iloadmax,maxLoadCondition,maxgamma,diagramIntAxis1,...
    rotCoordDiscISR,rotSection,cmax,CP]=EffISRColsRot(b,h,t,load_conditions,...
    rec,fy,npdiag,fdpc,E,beta1) 
%------------------------------------------------------------------------
% Syntax:
% [maxef,tablaEff,iloadmax,maxLoadCondition,maxgamma,diagramIntAxis1,...
%  rotCoordDiscISR,rotSection,cmax,CP]=EffISRColsRot(b,h,t,...
%  load_conditions,rec,fy,npdiag,fdpc,E,beta1)
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
% 
%------------------------------------------------------------------------
% PURPOSE: To determine the structural efficiency of a rectangular column
% cross-section reinforced with an Idealized Steel Profile, for each of the
% given load combinations.
%
% Note: the cross-section is rotated according to each load combinations
% of biaxial bending moments so that an interaction diagram is computed for
% each of them, respectively, with respect to their resultant bending
% moment axis of action.
% 
% OUTPUT: maxef:               Critical structural efficiency
%
%         tablaEff:            Resume table of the rebar design's 
%                              structural efficiency for each given load 
%                              condition
%
%         iloadmax:            is the index of the critical load condition
%                              according to the order in which they were
%                              given
%
%         cmax:                is the neutral axis depth of the critical
%                              load condition
%
%         CP:                  Rotated cross-section's plastic center 
%                              depth (see Doc) 
%
%         maxLoadCondition:    is the critical load condition in format:
%                              [n-load, Pu, Mu], where:
%                              Mu = sqrt( Mx^2 + My^2 )
%
%         diagramIntAxis1:     interaction diagram's data. Format:
%                              [ Pu, Mu, Mr, Fr * Mr, ecc ]
%
%         rotCoordDiscISR:     Rotated discrete ISR's elements' coordinates
%                              with respect to the most critical load 
%                              condition 
%
%         rotSection:          Are the coordinates for rotated
%                              cross-section's corners, with respect to the
%                              most critical load condition
% 
% INPUT:  b,h:                cross-section dimensions of column (width 
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
%         E:                  Elasticity modulus of the reinforcing steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

maxef=-inf;
nloads=length(load_conditions(:,1));
for i=1:nloads
    %% Interaction diagram in the load direction
    Pu=load_conditions(i,2);
    Mux=load_conditions(i,3);
    Muy=load_conditions(i,4);

    [diagramIntAxis,pot,poc,cvector1,newdispositionRebar,...
    newCoordCorners,newCP,gamma]=InteracDiagAxisISR...
    (npdiag,t,beta1,[b,h],rec,fy,fdpc,E,Mux,Muy);

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
        rotCoordDiscISR=newdispositionRebar;
        rotSection=newCoordCorners;
        cmax=c;
        CP=newCP;
        maxLoadCondition=trans_load_condition;
    end
end