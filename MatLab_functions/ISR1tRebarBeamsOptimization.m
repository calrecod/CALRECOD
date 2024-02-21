function [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,arrangement_t1,...
    arrangement_t2,maxEf,bestMr,area_var_t]=ISR1tRebarBeamsOptimization(E,b,...
    h,fy,fc,b_rec,h_rec,conditions,t2,pu_beams,wac,span,varAvailable)
    
%------------------------------------------------------------------------
% Syntax:
% [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,arrangement_t1,...
%    arrangement_t2,maxEf,bestMr,area_var_t]=ISR1tRebarBeamsOptimization(E,b,...
%    h,fy,fc,b_rec,h_rec,conditions,t2,pu_beams,wac,span,varAvailable)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To design optimally a rebar distribution over a beam 
% cross-section given the ISR.
% 
% OUTPUT: sepbarsRestric:       Is the parameter that indicates if the 
%                               minimum rebar separation restriction of 
%                               rebars in tension for the cross-section in 
%                               question is being complied: (1) indicates 
%                               that such restriction is not being complied,
%                               (0) indicates that such restriction is being
%                               complied
%
%         cbest:                is the depth of the neutral axis for the 
%                               optimized beam reinforced cross-section
%                               considering the optimal rebar design
%
%         b,h:                  are the final cross-section dimensions in
%                               case they suffered modifications after the
%                               optimal design process
%
%         bestBarDisposition:   are the local coordinates of rebar 
%                               disposition over the optimal designed 
%                               cross-section 
%
%         arrangement_t1:       are the list of rebar type transformed from
%                               the ISR in tension: a vector consisting of
%                               one column of length nbars in tension
%
%         arrangement_t2:       are the list of rebar type transformed from 
%                               the ISR in compression: a vector consisting 
%                               of one column of length nbars in compression
%
%         maxEf:                is the optimal final structural efficiency 
%                               for the optimal designed beam cross-section 
%                               considering the optimal rebar
%
%         bestMr:               is the optimal final bending resistance for
%                               the optimal designed beam cross-section 
%                               considering the optimal rebar
%
%         area_var_t:           is a vector consisting of the total optimal 
%                               rebar area in tension and compression (it 
%                               can be considered later for the assessment 
%                               of modification of inertia as a cracked 
%                               section)
%
% INPUT:  b_rec,h_rec:          are the concrete cover parameters horizontally and
%                               vertically, respectively
%
%         fc:                   is the f'c
%
%         E:                    is the Elasticity Modulus of rebar 
%
%         fy:                   is the yield stress of reinforcing steel
%
%         t2:                   is the optimal ISR consisting of vector of
%                               two elements [t_tension,t_compression]
%
%         pu_beams:             is the unitary cost of rebar assembly in 
%                               beams considering an average of assembly 
%                               performance (it is considered that various 
%                               types of rebars are placed simultaneously
%                               in the beam element along its length)
%
%         conditions:           Is the array containing the load conditions
%                               in format: [nload, Mu]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

if fc<2000 % units: Kg,cm
    betac=1.05-fc/1400;
    if (betac<0.65)
        betac=0.65;
    elseif (betac>0.85)
        betac=0.85;
    end
else        % units: lb, in
    betac=0.85-0.05*(fc-4000)/1000;
    if betac<0.65
        betac=0.65;
    elseif betac>0.85
        betac=0.85;
    end
end
fdpc=0.85*fc;

% 1. Determine all the possible reinforcement combinations _____________ 
[sepbarsRestric,cbest,b,h,maxEf,bestMr,area_var_t,nv_t,arrangement_t1,arrangement_t2,list_pac_t1,...
    list_pac_t2,bestBarDisposition]=RebarOptimalDesignBeams(b,h,b_rec,h_rec,...
    varAvailable,t2,fdpc,E,fy,conditions,betac);       

bestCost = EvaluateCostbeams(nv_t(1),nv_t(2),arrangement_t1,...
            arrangement_t2,pu_beams,varAvailable,wac,span);
