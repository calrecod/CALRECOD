function [maxef,tabEff,cxy]=effRecColsLinearSearch(diagrama,...
                            load_conditions,pot,poc,c_vector_bar)
%------------------------------------------------------------------------
% Syntax:
% [maxef,tabEff,cxy]=effRecColsLinearSearch(diagrama,...
%                             load_conditions,pot,poc,c_vector_bar)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the structural efficiency of a symmetrically 
%          reinforced column cross-section, given multiple load conditions.
%          The Bresler's formula is used (for biaxial bending compression)
%          and the Euclidean distance is used (for uniaxial bending
%          compression).
% 
% OUTPUT: maxef:                is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         tabEff:               is a table containing the structural 
%                               efficiency analysis data: size = [nload,8],
%                               in format: 
%           _____________________________________________________
%           [Pu,Mux,Muy,P{Rx},M{Rx},P{Ry},M{Ry},efficiency]
%           _____________________________________________________
%
%         cxy:                  are the neutral axis depth values 
%                               corresponding to the critical load condition, 
%                               for both axis directions: [cx,cy]
%
% INPUT:  diagrama:             is the interaction diagram data
%
%         load_conditions:      is the array containing the load conditions:
%                               size = [nload,4] in format [nload,Pu,Mux,Muy]
%
%         pot,poc:              are the max resistant axial force in tension
%                               of reinforcement steel (concrete is not 
%                               considered) and compression of the whole 
%                               reinforced cross-section area (both concrete
%                               area and rebar area are considered) 
%
%         c_vector_bar:         is the array containing the neutral axis 
%                               depth values for both cross-section axis 
%                               directions for all interaction diagram 
%                               points: size = (npoints+2) x 2
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

%% RESISTANCE COMPUTATION 
cx=[];
cy=[];
nconditions=length(load_conditions(:,1));
for sentido=1:2
    for j=1:nconditions
        pu=load_conditions(j,2);
        mu=load_conditions(j,2+sentido);
        penu=1/(mu/pu); % NOTE: The slope is the one taken for comparison
        
        tabEff(j,1)=pu;
        tabEff(j,1+sentido)=mu;
        
        k=1;

        ydos=diagrama(k,4*sentido-1);
        ytres=diagrama(k+1,4*sentido-1);
        xdos=diagrama(k,4*sentido);
        xtres=diagrama(k+1,4*sentido);
        pendos=ydos/xdos;
        pentres=ytres/xtres;
        while pendos<penu && pentres<penu
            k=k+1;
            ydos=diagrama(k,4*sentido-1);
            ytres=diagrama(k+1,4*sentido-1);
            xdos=diagrama(k,4*sentido);
            xtres=diagrama(k+1,4*sentido);
            pendos=ydos/xdos;
            pentres=ytres/xtres;
            if pendos<penu && pentres>penu
                break;
            end
        end
        
        mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
            (pu/mu-((ytres-ydos)/(xtres-xdos)));
        if mr>max(diagrama(:,4*sentido))
            [mr,k]=max(diagrama(:,4*sentido));
        end
        cdos=c_vector_bar(k,sentido);
        ctres=c_vector_bar(k+1,sentido);
        
        c=(cdos+ctres)*0.5;
        if sentido==1
            cx=[cx;c];
        elseif sentido==2
            cy=[cy;c];
        end
        if mr==0
            mr=0.00001;
        end
        pr=pu/mu*mr;

        tabEff(j,2+2*sentido)=pr;
        tabEff(j,3+2*sentido)=mr;
    end
end

%% EFFICIENCY COMPUTATION
maxef=-inf;
imax=1;
for j=1:nconditions
    prx=tabEff(j,4);
    pry=tabEff(j,6);
    mrx=tabEff(j,5);
    mry=tabEff(j,7);
    mux=tabEff(j,2);
    muy=tabEff(j,3);
    
    % Breler's formula
    % Inverse load equation
    prR=1/(1/abs(prx)+1/abs(pry)-1/pot);
    efpr=abs(prR/poc);
    % Countour load equation
    efmr=(abs(mux/mrx))^0.8+(abs(muy/mry))^0.8;
    if (efpr<0.1)
        ef=efmr;
    else
        ef=efpr;
    end
    
    % When only subject to uniaxial bending-compression then
    % the euclidean distance formula is used for the 
    % computation of the resistance efficiency
    if abs(muy)<0.1
        ef=sqrt(pu^2+mux^2)/sqrt(prx^2+mrx^2);
    elseif abs(mux)<0.1
        ef=sqrt(pu^2+muy^2)/sqrt(pry^2+mry^2);
    end
    if (ef>=maxef)
        imax=j;
        maxef=ef;
    end
    tabEff(j,8)=ef;
end
% To store the neutral axis depth values corresponding to the critical
% load condition
cxy=[cx(imax) cy(imax)];
