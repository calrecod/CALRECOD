function [maxef,diagramaInteraccion,tablaEficiencias,cxy]=eficienciaRecISRCols...
    (nconditions,diagramaInteraccion,load_conditions,pot,poc,c_vector)

%------------------------------------------------------------------------
% Syntax:
% [maxef,diagramaInteraccion,tablaEficiencias,cxy]=eficienciaRecISRCols...
%    (nconditions,diagramaInteraccion,load_conditions,pot,poc,c_vector)
%
%------------------------------------------------------------------------
% PURPOSE: The function computes the interaction diagram of a column 
% cross-section idealized with an 1t-ISR and its structural efficiency 
% based on the given load conditions.
% 
% OUTPUT: maxef:               is the critical structural efficiency of the
%                              column cross-section given different load 
%                              conditions
%
%         diagramaInteraccion: is the interaction diagram for the reinforced
%                              column cross-section
%
%         tablaEficiencias:    is the resume table of results consisting of
%                              nloads rows and eight columns as: 
%                              [Pu,Mux,Muy,PRx,PRy,MRx,MRy,Eff]
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
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

ncondiciones=length(load_conditions(:,1));
npuntos=length(diagramaInteraccion(:,1));
tablaEficiencias=zeros(ncondiciones,8);

for sentido=1:2
    if sentido==1
        mrxmax=max(diagramaInteraccion(:,4*sentido));
    elseif sentido==2
        mrymax=max(diagramaInteraccion(:,4*sentido));
    end
    % RESISTANCE COMPUTATION ___________________________________

    for j=1:nconditions
        pu=load_conditions(j,2);
        mu=load_conditions(j,2+sentido);
        k=1;
        tablaEficiencias(j,1)=pu;
        tablaEficiencias(j,1+sentido)=mu;
        mcondicion=mu/pu;
        
        ydos=diagramaInteraccion(k,4*sentido-1);
        xdos=diagramaInteraccion(k,4*sentido);
        ytres=diagramaInteraccion(k+1,4*sentido-1);
        xtres=diagramaInteraccion(k+1,4*sentido);
        
        cdos=c_vector(k);
        ctres=c_vector(k+1);
        mtres=diagramaInteraccion(k+1,4*sentido+1);
        if(mcondicion>0)
            mtres=diagramaInteraccion(k+1,4*sentido+1);
            while(mtres<0)
                ydos=diagramaInteraccion(k,4*sentido-1);
                xdos=diagramaInteraccion(k,4*sentido);
                
                ytres=diagramaInteraccion(k+1,4*sentido-1);
                xtres=diagramaInteraccion(k+1,4*sentido);
                
                cdos=c_vector(k);
                ctres=c_vector(k+1);
                
                mdos=diagramaInteraccion(k,4*sentido+1);
                mtres=diagramaInteraccion(k+1,4*sentido+1);

                k=k+1;
            end

            while(mtres>mcondicion)
                ytres=diagramaInteraccion(k+1,4*sentido-1);
                xtres=diagramaInteraccion(k+1,4*sentido);

                ydos=diagramaInteraccion(k,4*sentido-1);
                xdos=diagramaInteraccion(k,4*sentido);
                mdos=diagramaInteraccion(k,4*sentido+1);
                mtres=diagramaInteraccion(k+1,4*sentido+1);  
                
                cdos=c_vector(k);
                ctres=c_vector(k+1);
                
                if(mdos>mcondicion && mtres<mcondicion)
                    break;
                end

                k=k+1;
                if k>npuntos-1
                    break;
                end
            end
        elseif mcondicion<0
            k=1;
            while(mtres>mcondicion)
                ytres=diagramaInteraccion(k+1,4*sentido-1);
                xtres=diagramaInteraccion(k+1,4*sentido);

                ydos=diagramaInteraccion(k,4*sentido-1);
                xdos=diagramaInteraccion(k,4*sentido);
                mdos=diagramaInteraccion(k,4*sentido+1);
                mtres=diagramaInteraccion(k+1,4*sentido+1);   
                
                cdos=c_vector(k);
                ctres=c_vector(k+1);
                
                if(mdos>mcondicion && mtres<mcondicion)
                    break;
                end
                k=k+1;
                if k>npuntos-1
                    break;
                end
            end
        end
       
        mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
            (pu/mu-((ytres-ydos)/(xtres-xdos)));
        
        c=(cdos+ctres)*0.5;
        
        if sentido==1
            if mr>mrxmax
                mr=mrxmax;
            end
            cx=c;
        elseif sentido==2
            if mr>mrymax
                mr=mrymax;
            end
            cy=c;
        end
        pr=pu/mu*mr;

        tablaEficiencias(j,2+2*sentido)=pr;
        tablaEficiencias(j,3+2*sentido)=mr;
    end
end

cxy=[cx cy];
% EFFICIENCY COMPUTATION__________________________________________
maxef=0;
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
        maxef=ef;
    end
    tablaEficiencias(j,8)=ef;
end
