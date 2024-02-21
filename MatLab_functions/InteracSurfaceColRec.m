function [supX,supY,supZ]=InteracSurfaceColRec(b,h,...
    comborebar,npdiag,fy,fdpc,beta1,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
    dispositionRebar)

%------------------------------------------------------------------------
% Syntax:
% [supX,supY,supZ]=InteracSurfaceColRec(b,h,...
%  comborebar,npuntos,fy,fdpc,beta1,E,numberRebars1,...
%  numberRebars2,numberRebars3,numberRebars4,rebarAvailable,...
%  dispositionRebar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: to compute the 3D interaction surface of a reinforced
% rectangular column cross-section with symmetrical or asymmetrical rebar.
% 
% OUTPUT: supX,supY,supZ:   are the points of the interaction surface in 
%                           3D system coordinate as (x,y,z) 
%
% INPUT:  fdpc:             is the factored value of the concrete
%                           compressive concrete f'c as 0.85 * f'c, 
%                           according to the ACI 318 code
%
%         b,h:              are the cross-section dimensions (width and
%                           height, respectively)
%
%         numberRebars1,
%         numberRebars2,
%         numberRebars3,
%         numberRebars4:    are the number of rebars to be placed for 
%                           each of the cross-section boundaries
%
%         beta1:            is determined as stablished by code, (see 
%                           Documentation)
%
%         dispositionRebar: are the local coordinates of rebars over 
%                           the cross-section
%
%         comborebar:       is the vector containing the four rebar
%                           diamemters' indices (from the given rebar 
%                           database table), in format:
%
%   [index-bar-upper, index-bar-lower, index-bar-left,index-bar-right]
%
%         npdiag:           number of points to be computed for the
%                           interaction diagrams at each angle step
%
%         fy:               is the yield stress of the reinforcing steel
%
%         E:                is the Modulus of Elasticity of the
%                           reinforcing steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

for i=1:359
    k=i+0.5;
    mx=cos(deg2rad(k));
    my=sin(deg2rad(k));
    
    [diagramIntAxis,pot,poc,cvector1,newdispositionRebar,...
    newCoordCorners,newCP]=InteractionDiagramAxis...
    (npdiag,comborebar,b,h,fy,fdpc,beta1,E,numberRebars1,...
    numberRebars2,numberRebars3,numberRebars4,...
    rebarAvailable,dispositionRebar,mx,my);
    
    supY(:,i)=diagramIntAxis(:,3);
    for j=1:npdiag
        supX(j,i)=abs(diagramIntAxis(j,4))*cos(deg2rad(k));
        supZ(j,i)=abs(diagramIntAxis(j,4))*sin(deg2rad(k));
    end
end

supX(:,360)=supX(:,1);
supY(:,360)=supY(:,1);
supZ(:,360)=supZ(:,1);

figure(2)
mesh(supX,supZ,supY)
xlabel('Mx')
ylabel('Mz')
zlabel('P')
title('Interaction surface of a rectangular column')
hold on
