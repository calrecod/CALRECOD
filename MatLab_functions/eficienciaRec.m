function [maxef,tablaEficiencias,cxy]=eficienciaRec(diagrama,...
                                    load_conditions,pot,poc,c_vector_bar)
                            
%------------------------------------------------------------------------
% Syntax:
% [maxef,tablaEficiencias,cxy]=eficienciaRec(diagramaInteraccion,...
%                             load_conditions,pot,poc,c_vector_bar)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the structural efficiency of a symmetrically 
%          reinforced column cross-section.
% 
% OUTPUT: maxef:                is the critical structural efficiency 
%                               corresponding to the critical load condition
%
%         tablaEficiencias:     is a table containing the structural 
%                               efficiency analysis data: size = [nload,8],
%                               in format: 
%           _____________________________________________________
%           [nload,Pu,Mux,Muy,P{Rx},M{Rx},P{Ry},M{Ry},efficiency]
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
%                               points: size = [npoints+2,2]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nconditions=length(load_conditions(:,1));
tablaEficiencias=zeros(nconditions,8);
npuntos=length(diagrama(:,1));
for sentido=1:2
    
    % RESISTANCE COMPUTATION ___________________________________

for j=1:nconditions
    pu=load_conditions(j,2);
    mu=load_conditions(j,2+sentido);
    k=1;
    tablaEficiencias(j,1)=pu;
    tablaEficiencias(j,1+sentido)=mu;
    mcondicion=mu/pu;
    
    mtres=diagrama(k+1,4*sentido+1);
    
    cdos=c_vector_bar(k,sentido);
    ctres=c_vector_bar(k+1,sentido);
        
    ydos=diagrama(k,4*sentido-1);
    ytres=diagrama(k+1,4*sentido-1);
    xdos=diagrama(k,4*sentido);
    xtres=diagrama(k+1,4*sentido);
    if(mcondicion>0)
        mtres=diagrama(k+1,4*sentido+1);
        while(mtres<0)
            ydos=diagrama(k,4*sentido-1);
            xdos=diagrama(k,4*sentido);
            ytres=diagrama(k+1,4*sentido-1);
            xtres=diagrama(k+1,4*sentido);

            cdos=c_vector_bar(k,sentido);
            ctres=c_vector_bar(k+1,sentido);
            
            mdos=diagrama(k,4*sentido+1);
            mtres=diagrama(k+1,4*sentido+1);
            
            k=k+1;
            if k>npuntos-1
                break;
            end
        end
    end
            
    while(mtres>mcondicion)
        ytres=diagrama(k+1,4*sentido-1);
        xtres=diagrama(k+1,4*sentido);
 
        ydos=diagrama(k,4*sentido-1);
        xdos=diagrama(k,4*sentido);
        mdos=diagrama(k,4*sentido+1);
        mtres=diagrama(k+1,4*sentido+1); 
        
        cdos=c_vector_bar(k,sentido);
        ctres=c_vector_bar(k+1,sentido);
            
        if(mdos>mcondicion && mtres<mcondicion)
            break;
        end

        k=k+1;
        if k>npuntos-1
            break;
        end
    end
       
    mr=(((ydos-ytres)/(xtres-xdos))*xtres+ytres)/...
        (pu/mu-((ytres-ydos)/(xtres-xdos)));
    
    c=(cdos+ctres)*0.5;
    if sentido==1
        cx=c;
    elseif sentido==2
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


