function [maxef,tablaEficiencias,cxy]=effRecColsBinarySearch...
        (diagramaInteraccion,load_conditions,pot,poc,c_vector)

%------------------------------------------------------------------------
% Syntax:
% [maxef,tablaEficiencias,cxy]=eficienciaRecColsBinarySearch...
%    (diagramaInteraccion,load_conditions,pot,poc,c_vector)
%
%------------------------------------------------------------------------
% PURPOSE: The function determines the structural efficiency of a rectangular
% column cross-section subject to multiple load biaxial combinations.
% Binary search is deployed to determine the corresponding resistance for
% every load combination based on its eccentricity.
% 
% OUTPUT: maxef:               is the critical structural efficiency of the
%                              column cross-section given different load 
%                              conditions
%
%         tablaEficiencias:    is the resume table of results consisting of
%                              nloads rows and eight columns as: 
%                              [Pu,Mux,Muy,PRx,MRx,PRy,MRy,Eff]
%
%         cxy:                 is a vector containing the neutral axis depth
%                              of each cross-section direction according to
%                              the most critical load condition as [cx,cy]
%
% INPUT:  poc,pot:             are the resistant axial load in compression 
%                              and axial load in tension, respectively, of
%                              the reinforced column cross-section
%
%         load_conditions:     is the vector consisting of nload_combinations
%                              as rows and four columns corresponding to 
%                              [nload,Pu,Mux,Muy]
%
%         diagramaInteraccion: is the array containing the points of the 
%                              interaction diagram for both axis directions
%                              of the column cross-section as (Pu, Muaxis)
%                              determined with function widthEfficiencyCols
%                              (see Documentation)
%
%         c_vector:            is an array consisting of n interaction 
%                              diagram points as rows and two columns. Each
%                              combination [cxi,cyi] corresponds to each 
%                              point of the interaction diagram which is 
%                              determined with the function bisectionMr4t 
%                              (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-07-19
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nconditions=length(load_conditions(:,1));
npuntos=length(diagramaInteraccion(:,1));
tablaEficiencias=zeros(nconditions,8);
CX=[];
CY=[];
for sentido=1:2
    % RESISTANCE COMPUTATION ___________________________________
    
    for j=1:nconditions
        pu=load_conditions(j,2);
        mu=load_conditions(j,2+sentido);
        penu=1/(mu/pu); % m=1/e NOTE: The slope is the one used for 
                        % comparison
        
        tablaEficiencias(j,1)=pu;
        tablaEficiencias(j,1+sentido)=mu;
        
        k=1;
        n1=1;
        nmiddle=round(npuntos/2);
        n2=npuntos;
        
        ytres=diagramaInteraccion(nmiddle,4*sentido-1);
        xtres=diagramaInteraccion(nmiddle,4*sentido);
        p3=ytres/xtres;
        dif=inf;
        if nmiddle==1
            nmiddle=2;
        end
        while dif>=2
            if penu>p3
                ydos=diagramaInteraccion(nmiddle-1,4*sentido-1);
                xdos=diagramaInteraccion(nmiddle-1,4*sentido);

                ytres=diagramaInteraccion(nmiddle,4*sentido-1);
                xtres=diagramaInteraccion(nmiddle,4*sentido);
                p3=ytres/xtres;
                cdos=c_vector(nmiddle-1,sentido);
                ctres=c_vector(nmiddle,sentido);
                dif=n2-n1;
                
                n1=nmiddle;
                n2=n2;
                nmiddle=nmiddle+round(1/(2*2^k)*npuntos);
                
            elseif penu<p3
                ydos=diagramaInteraccion(nmiddle,4*sentido-1);
                xdos=diagramaInteraccion(nmiddle,4*sentido);

                ytres=diagramaInteraccion(nmiddle+1,4*sentido-1);
                xtres=diagramaInteraccion(nmiddle+1,4*sentido);
                p3=ytres/xtres;
                cdos=c_vector(nmiddle,sentido);
                ctres=c_vector(nmiddle+1,sentido);
                dif=n2-n1;
                
                n1=n1;
                n2=nmiddle;
                nmiddle=nmiddle-round(1/(2*2^k)*npuntos);
                
            end
            k=k+1;
        end
        mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
            (pu/mu-((ytres-ydos)/(xtres-xdos)));
        c=(cdos+ctres)*0.5;
        if sentido==1
            cx=c;
            CX=[CX;
                cx];
        elseif sentido==2
            cy=c;
            CY=[CY;
                cy];
        end
        pr=pu/mu*mr;

        tablaEficiencias(j,2+2*sentido)=pr;
        tablaEficiencias(j,3+2*sentido)=mr;
    end
end

% EFFICIENCY COMPUTATION__________________________________________
maxef=0;
imax=0;
for j=1:nconditions
    prx=tablaEficiencias(j,4);
    pry=tablaEficiencias(j,6);
    mrx=tablaEficiencias(j,5);
    mry=tablaEficiencias(j,7);
    mux=tablaEficiencias(j,2);
    muy=tablaEficiencias(j,3);
    
    prR=1/(1/prx+1/pry-1/pot);
    efpr=prR/poc;
    efmr=(mux/mrx)+(muy/mry);
    if (efpr<0.1)
        ef=efmr;
    else
        ef=efpr;
    end
    if (maxef<ef)
        imax=j;
        maxef=ef;
    end
    tablaEficiencias(j,8)=ef;
end
cxy=[CX(imax) CY(imax)];