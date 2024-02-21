function [ef,tabEff,c]=effColsRot1DirecLS(diagrama,...
                           resultLoadConditions,c_vector_bar)
%------------------------------------------------------------------------
% Syntax:
% [ef,tabEff,c]=effColsRot1DirecLS(diagrama,...
%                   resultLoadConditions,c_vector_bar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%------------------------------------------------------------------------
% PURPOSE: To compute the structural efficiency of a rotated reinforced
%          column of rectangular cross-section, given multiple load
%          conditions. The function deploys linear search to find the
%          corresponding resistance for multiple load combinations 
%          according to their load eccentricity from the interaction 
%          diagram's data.
% 
% OUTPUT: ef:                   is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         tablaEficiencias:     is a table containing the structural 
%                               efficiency analysis data: size = [nload,8],
%                               in format: 
%                   _______________________________
%                       [Pu,Mu,PR,MR,efficiency]
%                   _______________________________
%
%         c:                    is the neutral axis depth value
%                               corresponding to the most critical load 
%                               condition for the rotated cross-section
%
% INPUT:  diagrama:             is the interaction diagram data
%                   _______________________________
%                       [Pn,Mn,PR,MR,eccentricity]
%                   _______________________________
%
%         resultLoadConditions: is the array containing the load conditions:
%                               size = [nload,4] in format [nload,Pu,Mux,Muy]
%
%         c_vector_bar:         is the array containing the neutral axis 
%                               depth values for both cross-section axis 
%                               directions for all interaction diagram 
%                               points: size = [npoints+2,2]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
%% RESISTANCE COMPUTATION
pu=resultLoadConditions(1,2);
mu=resultLoadConditions(1,3);
penu=1/(mu/pu); % NOTE: The slope is the one taken for comparison

tabEff(1,1)=pu;
tabEff(1,2)=mu;

k=1;
npdiag=length(diagrama(:,1));
ydos=diagrama(k,3);
ytres=diagrama(k+1,3);
xdos=diagrama(k,4);
xtres=diagrama(k+1,4);
pendos=ydos/xdos;
pentres=ytres/xtres;
while pendos<penu && pentres<penu
    k=k+1;
    ydos=diagrama(k,3);
    ytres=diagrama(k+1,3);
    xdos=diagrama(k,4);
    xtres=diagrama(k+1,4);
    pendos=ydos/xdos;
    pentres=ytres/xtres;
    if pendos<penu && pentres>penu
        break;
    end
    if k==npdiag-1
        break;
    end
end

mr=abs((((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
    (pu/mu-((ytres-ydos)/(xtres-xdos))));

if mr==0 || isnan(mr)
    mr=1e-6; % a small value of mr in case mu=0 but not equal...
             % to zero to avoid division by 0
elseif mr>max(diagrama(:,4))
    [mr,k]=max(diagrama(:,4));
end
pr=pu/mu*mr;

cdos=c_vector_bar(k,1);
ctres=c_vector_bar(k+1,1);

c=(cdos+ctres)*0.5;

tabEff(1,3)=pr;
tabEff(1,4)=mr;

%% COMPUTATION OF STRUCTURAL EFFICIENCY
ef=sqrt(pu^2+mu^2)/sqrt(pr^2+mr^2);
if isnan(ef)==1
    ef=1e10;
end
tabEff(1,5)=ef;
    