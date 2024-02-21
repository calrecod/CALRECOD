function [maxef,tablaEficiencias,c]=effCircColsLS(diagrama,...
                                    load_conditions,c_vector)
                            
%------------------------------------------------------------------------
% Syntax:
% [maxef,tablaEficiencias,c]=effCircColsLS(diagrama,...
%                                   load_conditions,c_vector)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the structural efficiency of a symmetrically 
%          reinforced column of circular cross-section, given multiple 
%          load conditions.
% 
% OUTPUT: maxef:                is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         tablaEficiencias:     is a table containing the structural 
%                               efficiency analysis data: size = [nload,8],
%                               in format: 
%                   _______________________________
%                       [Pu,Mu,PR,MR,efficiency]
%                   _______________________________
%
%         c:                    is the neutral axis depth values
%                               corresponding to the critical load condition
%
% INPUT:  diagrama:             is the interaction diagram data
%
%         load_conditions:      is the array containing the load conditions:
%                               size = [nload,4] in format [nload,Pu,Mux,Muy]
%
%         c_vector:             is the array containing the neutral axis 
%                               depth values for all interaction diagram 
%                               points: vector size: [npoints+2,2]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

nconditions=length(load_conditions(:,1));
tablaEficiencias=zeros(nconditions,5);
C=[];
    
% RESISTANCE COMPUTATION ___________________________________
for j=1:nconditions
    pu=load_conditions(j,2);
    mu=load_conditions(j,3);
    penu=1/(mu/pu); % NOTE: The slope is the one taken for comparison

    tablaEficiencias(j,1)=pu;
    tablaEficiencias(j,2)=mu;

    k=1;

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
    end
    mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
        (pu/mu-((ytres-ydos)/(xtres-xdos)));
    if mr>max(diagrama(:,4))
        [mr,k]=max(diagrama(:,4));
    end
    cdos=c_vector(k,1);
    ctres=c_vector(k+1,1);
    
    c=(cdos+ctres)*0.5;
    
    C=[C;
       c];
    pr=pu/mu*mr;

    tablaEficiencias(j,3)=pr;
    tablaEficiencias(j,4)=mr;

end
%% EFFICIENCY COMPUTATION
maxef=-10000;
imax=0;
for j=1:nconditions
    pr=tablaEficiencias(j,3);
    mr=tablaEficiencias(j,4);
    mu=tablaEficiencias(j,2);
    
    ef=sqrt(pu^2+mu^2)/sqrt(pr^2+mr^2);
    tablaEficiencias(j,5)=ef;
    
    if (ef>maxef)
        imax=j;
        maxef=ef;
    end
    
end
c=C(imax);