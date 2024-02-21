function [ISReMec]=eleMecISRSymRecRot(discreteISRcoord,ndt,da4,b,h,c,...
         fdpc,E,beta1,CP,RotCornerSec,gamma)
%------------------------------------------------------------------------
% Syntax:
% [ISReMec]=eleMecISRSymRecRot(discreteISRcoord,ndt,da4,b,h,c,...
%  fdpc,E,beta1,CP,RotCornerSec,gamma)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the sum of resistant forces of a rotated 
% asymmetrically reinforced column cross-section considering an Idealized
% Smeared Reinforcement ISR and the concrete zone in compression.
% 
% OUTPUT: ISReMec:              vector that contains the output [Fs,Ms;
%                                                                Fc,Mc]
%
% INPUT:  E:                    Modulus of Elasticity of steel
%                               reinforcement
%
%         b,h:                  cross-section dimensions
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta1:                is determined as stablished by code 
%                               (see Documentation)
%
%         ndt:                  is the total number of discrete ISR
%                               elements over the cross-section. 
%
%         discreteISRcoord:     are the local coordinates of ISR's discrete
%                               element over the cross-section
%
%         CP:                   Plastic Center location for each axis
%                               direction of the column cross-section
%
%         RotCornerSec:         coordinates of each of the four
%                               rotated cross-section's corners
%
%         gamma:                angle of rotation for the cross-section
%
%         da4:                  transversal cross-section area of each
%                               discrete ISR element at each of the four
%                               cross-section's boundaries
%
%         c:                    is the neutral axis's depth
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

coordyCP=max(RotCornerSec(:,2))-CP;

eMecISR=zeros(ndt,9);
eMecISR(:,2)=discreteISRcoord(:,1);
eMecISR(:,3)=discreteISRcoord(:,2);

sumaM=0;
sumaF=0;
  
av1=da4(1);
av2=da4(2);
av3=da4(3);
av4=da4(4);

ndtb=ndt/4;
for i=1:ndt
    if i<=ndtb
        av=av1;
        eMecISR(i,5)=av;
    elseif i<=(2*ndtb) && i>(ndtb)
        av=av2;
        eMecISR(i,5)=av;
        
    elseif i<=(3*ndtb) && i>(2*ndtb)
        av=av3;
        eMecISR(i,5)=av;
        
     elseif i<=(4*ndtb) && i>(3*ndtb)
        av=av4;
        eMecISR(i,5)=av;   
        
    end
    % to calculate eps:
    eMecISR(i,4)=max(RotCornerSec(:,1))-eMecISR(i,3); 
    
    eMecISR(i,6)=0.003/c*(eMecISR(i,4)-c);
    if (eMecISR(i,6)<-0.0021)
        eMecISR(i,6)=-0.0021;
    elseif(eMecISR(i,6)>0.0021)
        eMecISR(i,6)=0.0021;
    end
    eMecISR(i,7)=eMecISR(i,6)*E;
    eMecISR(i,8)=eMecISR(i,7)*av;
    eMecISR(i,9)=eMecISR(i,8)*(coordyCP-eMecISR(i,3));
    
    sumaF=sumaF+eMecISR(i,8);
    sumaM=sumaM+eMecISR(i,9);
end

a=beta1*c;
elemConc=casoConcretoRecRot(a,fdpc,b,h,RotCornerSec,CP,gamma);
ISReMec=[sumaF sumaM;
    elemConc(1) elemConc(2)];