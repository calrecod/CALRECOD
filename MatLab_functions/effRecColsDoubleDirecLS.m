function [maxef,tableEff,cxy,Mr]=effRecColsDoubleDirecLS(diagrama1,...
           diagrama2,load_conditions,pot,poc,c_vector_bar1,c_vector_bar2)

%------------------------------------------------------------------------
% Syntax:
% [maxef,tableEff,cxy]=effRecColsDoubleDirecLS(diagrama1,...
%  diagrama2,load_conditions,pot,poc,c_vector_bar,c_vector_bar2)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the structural resistance efficiency of an 
%          asymmetrically reinforced rectangular column's cross-section
%          given certain load conditions. The function accepts both
%          positive and negative bending moments over a single axis
%          direction. A Linear-Search method is used to look for the
%          corresponding resistance of each load condition according to its
%          respective eccentricity
% 
% OUTPUT: maxef:                is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         tableEff:             is a table containing the structural 
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
% INPUT:  diagrama1:            is the first interaction diagram's data
%                               (positive bending moments)
%
%         diagrama2:            is the second interaction diagram's data
%                               (negative bending moments)
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
%         c_vector_bar1,
%         c_vector_bar2:        is the array containing the neutral axis 
%                               depth values for both cross-section axis 
%                               directions for all interaction diagram 
%                               points: size = [npoints,2] (x,y). Such
%                               vectors are obtained from the function:
%                               "EvalAsymDoubleDirection"
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-08-02
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% RESISTANCE COMPUTATION 
npdiag=length(diagrama1(:,1));
nconditions=length(load_conditions(:,1));
tableEff=zeros(nconditions,8);
CX=[];
CY=[];
for sentido=1:2
    for j=1:nconditions
        pu=load_conditions(j,2);
        mu=load_conditions(j,2+sentido);
        penu=1/(mu/pu); % NOTE: The slope is the one taken for comparison
        
        tableEff(j,1)=pu;
        tableEff(j,1+sentido)=mu;
        
        if mu>=0
            k=1;
            ydos=diagrama1(k,4*sentido-1);
            ytres=diagrama1(k+1,4*sentido-1);
            xdos=diagrama1(k,4*sentido);
            xtres=diagrama1(k+1,4*sentido);
            pendos=ydos/xdos;
            pentres=ytres/xtres;
            while pendos<penu && pentres<penu
                k=k+1;
                ydos=diagrama1(k,4*sentido-1);
                ytres=diagrama1(k+1,4*sentido-1);
                xdos=diagrama1(k,4*sentido);
                xtres=diagrama1(k+1,4*sentido);
                pendos=ydos/xdos;
                pentres=ytres/xtres;
                if pendos<penu && pentres>penu
                    break;
                end
                if k==npdiag-1
                    break;
                end
            end
            cdos=c_vector_bar1(k,sentido);
            ctres=c_vector_bar1(k+1,sentido);
            mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
                (pu/mu-((ytres-ydos)/(xtres-xdos)));
        elseif mu<0
            k=1;
            ydos=diagrama2(k,4*sentido-1);
            ytres=diagrama2(k+1,4*sentido-1);
            xdos=diagrama2(k,4*sentido);
            xtres=diagrama2(k+1,4*sentido);
            pendos=ydos/xdos;
            pentres=ytres/xtres;
            while pendos>penu && pentres>penu
                k=k+1;
                ydos=diagrama2(k,4*sentido-1);
                ytres=diagrama2(k+1,4*sentido-1);
                xdos=diagrama2(k,4*sentido);
                xtres=diagrama2(k+1,4*sentido);
                pendos=ydos/xdos;
                pentres=ytres/xtres;
                if pendos>penu && pentres<penu
                    break;
                end
                if k==npdiag-1
                    break;
                end
            end
            cdos=c_vector_bar2(k,sentido);
            ctres=c_vector_bar2(k+1,sentido);
            mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
                (pu/mu-((ytres-ydos)/(xtres-xdos)));
        end
        c=(cdos+ctres)*0.5;
        if sentido==1
            CX=[CX;
                c];
        elseif sentido==2
            CY=[CY;
                c];
        end
        pr=pu/mu*mr;

        tableEff(j,2+2*sentido)=pr;
        tableEff(j,3+2*sentido)=mr;
        
    end
end

%% STRUCTURAL EFFICIENCY COMPUTATION
maxef=0;
for j=1:nconditions
    prx=tableEff(j,4);
    pry=tableEff(j,6);
    mrx=tableEff(j,5);
    mry=tableEff(j,7);
    mux=tableEff(j,2);
    muy=tableEff(j,3);
    
    % Breler's formula
    % Inverse load equation
    prR=1/(1/abs(prx)+1/abs(pry)-1/abs(pot));
    efpr=prR/poc;
    % Countour load equation
    efmr=abs(mux/mrx)+abs(muy/mry);
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
    
    if (ef>maxef)
        imax=j;
        maxef=ef;
    end
    tableEff(j,8)=ef;
end
% To store the neutral axis depth values corresponding to the critical
% load condition
cxy=[CX(imax) CY(imax)];
Mr=[tableEff(imax,5),tableEff(imax,7)];