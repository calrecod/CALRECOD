function [Inertia_xy_modif,b,h,bestArrangement,best_disposicion,cost_elem_col,...
    Ac_sec_elem,Ef_sec_col,Mr_col]=isrColumnsSymAsym(pu_cols,height,b,h,rec,fy,...
    fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,optimaPlot,...
    plotsISRdiagrams,plotRebarDesign)

%------------------------------------------------------------------------
% Syntax:
% [Inertia_xy_modif,b,h,bestArrangement,best_disposicion,cost_elem_col,...
%   Ac_sec_elem,Ef_sec_col,Mr_col]=isrColumnsSymAsym(pu_cols,height,b,h,rec,fy,...
%   fc,load_conditions,cols_sym_asym_isr,condition_cracking,ductility,optimaPlot,...
%   plotsISRdiagrams,plotRebarDesign)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement design, either with a pure
%          ISR or with symmetric rebar.
% 
% OUTPUT: Inertia_xy_modif:         momentum of inertia of the optimal 
%                                   reinforced cross-section for both axis
%                                   directions as [Ix,Iy] considering the
%                                   reinforcement with cracked or non-cracked
%                                   section mechanisms
%
%         b,h:                      are the final cross-section dimensions
%                                   in case of a need of modification to 
%                                   comply with the restrictions criteria
%
%         bestArrangement:          is the list of rebar types for each rebar
%                                   of the optima reinforcement option: 
%                                   size = [nbars,1] consisting of a number
%                                   from 1 to 7 by default
%
%         best_disposicion:         is the array consisting of the local
%                                   rebar coordinates over the cross-section
%                                   for the optimal rebar option
%
%         cost_elem_col:            is the cost of the optimal rebar option 
%                                   (only steel is considered)
%
%         Ac_sec_elem:              is the optimal rebar area
%
%         Ef_sec_col:               is the critical structural efficiency
%                                   for the optimal reinforcement option
%
%         Mr_col:                   are the resisting moments for both axis 
%                                   directions of the cross-section as
%                                   [M{Rx}, M{Ry}]
%
% INPUT:  b,h:                      initial given cross-section dimensions
%
%         pu_cols:                  is the database of reinforcement assembly
%                                   and construction unit cost: format by 
%                                   default: 
%    -----------------------------------------------------------------
%    pu_col=[PU{#4}, PU{#5}, PU{#6}, PU{#8}, PU{#9}, PU{#10}, PU{#12}]
%    -----------------------------------------------------------------
%
%         height:                   is the total length of the column
%
%         f'c:                      compressive concrete strength
%
%         fy:                       yield strength of reinforcement bars
%
%         load_conditions:          load condition array: size = [nloads,4],
%                                   in format: [nload,Pu,Mux,Muy]
%
%         ductility:                is demand ductility parameters, options
%                                   are (1) for low ductility section 
%                                   requirements,(2) for medium ductility, 
%                                   (3) for high ductility
%
%         rec:                      are the concrete cover values for both
%                                   cross-section axis directions: 
%                                   [coverX,coverY]
%
%         cols_sym_asym_isr:        is the reinforcement option parameters,
%                                   options are: ''Symmetric'' or ''ISR''
%
%         condition_cracking:       is the cracking mechanisms to be 
%                                   considered, options are: ''Cracked'' or
%                                   ''Non-cracked''
%
%         optimPlot:                is the parameters that indicates if the 
%                                   optima rebar area convergence is required
%                                   or not. Options are: (1) they are required,
%                                   (2) they are not required
%
%         plotsISRdiagrams:         is the parameters that indicates if the
%                                   optima ISR interaction diagrams are 
%                                   required or not. Options are: (1) they
%                                   are required, (2) they are not required
%
%         plotRebarResults:         Is the parameters that indicates if the
%                                   rebar design results are required or not.
%                                   Options are: (1) they are required, 
%                                   (2) they are not required
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

pu=load_conditions(1,2);
mux=load_conditions(1,3);
muy=load_conditions(1,4);

excentricity_x=abs(mux/pu)*100; %cm
excentricity_y=abs(muy/pu)*100; %cm

excentricity_xy=[excentricity_x,excentricity_y];

E=2e6;
npuntos=20;
fdpc=fc*0.85;

if fc<280
    betac=0.85;
elseif fc>=280
    betac=1.05-fc/1400;
    if (betac<0.65)
        betac=0.65;
    elseif (betac>0.85)
        betac=0.85;
    end
end
[b,h,cost_elem_col,Ac_sec_elem,Ef_sec_col,Mr_col,t_value_01,t_value_03,cxy]=...
    isrColumns(pu_cols,height,b,h,rec,fy,fc,load_conditions,ductility,optimaPlot,...
    plotsISRdiagrams);

if cols_sym_asym_isr=="Symmetric"
    % calculo de inercia como secci?n no-agrietada
    [Mr_col,h,Inertia_xy_modif,Ac_sec_elem,cost_elem_col,ov,nv,Ef_sec_col,...
         bestArrangement,best_disposicion]=optimalrebarColsSym(b,h,rec,...
         Ac_sec_elem,E,npuntos,fdpc,betac,pu_cols,load_conditions,condition_cracking,...
         plotRebarDesign);
    
    cost_elem_col=cost_elem_col*0.01*height;
    
elseif cols_sym_asym_isr=="ISR"
    [Inertia_xy_modif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,t_value_01,excentricity_xy,...
                            t_value_03,pu,cxy,condition_cracking,E);
    bestArrangement=[];
    best_disposicion=[];
elseif cols_sym_asym_isr=="Asymmetric"
    rebarOptimConv=0; % this parameters indicates if the optima
                      % rebar convergence is required or not. See
                      % documentation
    
    t4=[t_value_01,t_value_01,t_value_03,t_value_03];
    [Mr_col,h,Inertia_xy_modif,Ac_sec_elem,cost_elem_col,nv,Ef_sec_col,bestArrangement,...
    best_disposicion]=PSOAsymmetricRebar(t4,b,h,rec,Ac_sec_elem,npuntos,fdpc,...
    betac,pu_cols,load_conditions,condition_cracking,plotRebarDesign,rebarOptimConv);
end