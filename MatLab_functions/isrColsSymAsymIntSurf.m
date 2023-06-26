function [b,h,bestArrangement,bestDisposition,costElemCol,AcSecElm,...
    EfsecCol,Mrcol,bestLoad]=isrColsSymAsymIntSurf(puCols,height,b,h,rec,...
    fy,fc,load_conditions,cols_sym_asym_isr,ductility,wac,RebarAvailable,...
    optimaPlot,plotsISRdiagrams,plotRebarDesign)
%------------------------------------------------------------------------
% Syntax:
% [b,h,bestArrangement,bestDisposition,costElemCol,Ac_sec_elem,...
%  Ef_sec_col,Mr_col]=isrColsSymAsymIntSurf(puCols,height,b,h,rec,fy,...
%  fc,load_conditions,cols_sym_asym_isr,ductility,wac,RebarAvailable,...
%  optimaPlot,plotsISRdiagrams,plotRebarDesign)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement design, either with a pure
%          ISR, with symmetrical or asymmetrical rebar.
% 
% Note: When either symmetrical or asymmetrical rebar is wanted, the 
% Inverse Load method (Bresler's formula) and the Contour Load method is 
% used for the design of columns.
%
% OUTPUT: b,h:                      are the final cross-section dimensions
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
%         puCols:                   is the database of reinforcement 
%                                   assembly and construction unit cost
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

E=fy/0.0021; % Modulus of Elasticity of reinforcing steel
npdiag=30; % number of points to compute for the interaction diagram
fdpc=fc*0.85; % reduced compressive strength

if fc<2000 % units: Kg,cm
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
else       % units: lb,in
    betac=0.85-0.05*(fc-4000)/1000;
    if betac<0.65
        betac=0.65;
    elseif betac>0.85
        betac=0.85;
    end
end
[b,h,costElemCol,AcSecElm,EfsecCol,Mrcol,t_value_01,t_value_03,bestcxy,...
bestLoad]=isrColsInterSurf(puCols,height,b,h,rec,fy,fc,load_conditions,...
wac,ductility,optimaPlot,plotsISRdiagrams);

if cols_sym_asym_isr=="Symmetric"
    [Mrcol,h,AcSecElm,costElemCol,bestdiagram,bestnv,EfsecCol,...
    bestArrangement,bestDisposition,nv4,bestcxy,bestLoad]=...
    supOptimRebarSymIntSurf(b,h,rec,AcSecElm,E,npdiag,fdpc,betac,...
    puCols,load_conditions,wac,height,ductility,RebarAvailable,...
    plotRebarDesign);
    
elseif cols_sym_asym_isr=="ISR"
    bestArrangement=[];
    bestDisposition=[];
elseif cols_sym_asym_isr=="Asymmetric"
    
    [h,AcSecElm,EfsecCol,bestdiagram,bestArrangement,bestDisposition,...
    Mrcol,bestcxy,bestCP,costElemCol,bestLoad]=AssemRebarSymAsymIntSurf(b,...
    h,rec,AcSecElm,E,npdiag,fdpc,betac,height,wac,puCols,...
    RebarAvailable,load_conditions,ductility,plotRebarDesign);
end