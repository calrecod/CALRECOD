function [sepbarsRestric,cbest,b,h,maxEf,bestMr,area_var_t,nv_t,arreglo_t1,...
    arreglo_t2,list_pac_t1,list_pac_t2,disposition_rebar]=RebarOptimalDesignBeams...
    (b,h,b_rec,h_rec,varDisponibles,t2,fdpc,E,fy,load_conditions,betac)
    
%------------------------------------------------------------------------
% Syntax:
% [sepbarsRestric,cbest,b,h,maxEf,bestMr,area_var_t,nv_t,arreglo_t1,...
% arreglo_t2,list_pac_t1,list_pac_t2,disposition_rebar]=RebarOptimalDesignBeams...
% (b,h,b_rec,h_rec,varDisponibles,t2,fdpc,E,fy,load_conditions,betac)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal rebar arrangement based on minimum area
% for a beam cross-section.
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
%         disposition_rebar:    are the local coordinates of rebar 
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
%         nv_t                  Vector consisting of the number of rebars in
%                               tension and compression, stated as 
%                               [nrebar_tension,nrebar_compression]
%
%         list_pac_t1,          Vectors that contain a value either 1 or 2, 
%         list_pac_t2:          which indicate the manner in which a rebar 
%                               is displayed over the cross-section. Number
%                               1 means that the rebar is laid out at the 
%                               outer horizontal layer. Number 2 means that 
%                               the rebar is laid out at the inner horizontal
%                               layer (see Documentation)
%
% INPUT:  b_rec,h_rec:          are the concrete cover parameters horizontally and
%                               vertically, respectively (cm)
%
%         fdpc:                 is the f'c reduced by the factor 0.85 
%                               according to code, in units (Kg/cm2)
%
%         E:                    is the Elasticity Modulus of steel in 
%                               Kg/cm2
%
%         fy:                   is the yielding stress in Kg/cm2
%
%         t2:                   is the optimal ISR consisting of vector of
%                               two elements [t_tension,t_compression]
%
%         load_conditions:      Is the array containing the load conditions
%                               in format: [nload, Mu]
%
%         varDisponibles:       Commercial available type of rebars database
%                               as indicated in function ISR1tRebarBeamsOptimization
%                               (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

tma=2/3; %(in)
d=h-h_rec;
bp=b-2*b_rec;

ntypes=length(varDisponibles(:,1)); % tipos diferentes de varillas

at1=t2(1)*bp;
atmin=t2(2)*bp;

at2=[at1,atmin];
% to determine number of rebars for each available option
area_var_t=zeros(1,2);
nv_t=zeros(1,2);
for i=1:2

    ac_bar_space=zeros(ntypes,3); %[areabar espacio]

    amin=inf;
    for k=1:ntypes
        j=0;
        tipovar=k;
        areavarj=varDisponibles(tipovar,2)^2*pi/4;
        nv=fix(at2(i)/areavarj)+1;
        if nv==1
            nv=2;
        end
        diamj=varDisponibles(tipovar,2);
        atotal=nv*areavarj;
        % Rebar separation:
        if fdpc<2000 % units: (Kg,cm)
            sepMin=max([3/2*tma*2.54,diamj]);
        else         % units: (lb,in)
            sepMin=max([3/2*tma,diamj]);
        end
        
        sep1=(b-2*b_rec-nv*diamj)/(nv-1);

        sep2=(b-2*b_rec-(fix(nv*0.5)+1)*diamj)/(fix(nv*0.5));

        ac_bar_space(k,1)=atotal; 
        ac_bar_space(k,2)=sep1; % in one rebar pack
        ac_bar_space(k,3)=nv;
        ac_bar_space(k,4)=sep2; % in two rebar packs

        if ac_bar_space(k,1)<amin
            amin=atotal;
            index_min=k;
            best_area_var=areavarj;
            if sep1<sepMin && sep2>sepMin

                pac=2;
                sep_bars=sep2;
            elseif sep1>=sepMin

                sep_bars=sep1;
                pac=1;
            elseif sep1<sepMin && sep2<sepMin
                pac=2;
                sep_bars=sep2;
            end
        end

    end
    if sep_bars<sepMin
        sepbarsRestric=1;
    else
        sepbarsRestric=0;
    end
    nv=ac_bar_space(index_min,3);

    arreglo_var=zeros(nv,1)+index_min;
    if pac==1
        list_pac=zeros(1,nv)+pac;
    elseif pac==2
        if mod(nv,2)~=0
            list_pac=zeros(1,fix(nv/2)+1)+pac;
            list_pac(fix(nv/2)+1)=1;
        else
            list_pac=zeros(1,fix(nv/2))+pac;
        end
    end
    area_var_t(i)=nv*best_area_var;
    nv_t(i)=nv;
    if i==1
        arreglo_t1=arreglo_var;
        list_pac_t1=list_pac;

    elseif i==2
        arreglo_t2=arreglo_var;
        list_pac_t2=list_pac;
    end
end
Mu=load_conditions(1,2);

[disposition_rebar]=RebarDisposition1tBeams(b,h,b_rec,h_rec,...
varDisponibles,nv_t,arreglo_t1,arreglo_t2,list_pac_t1,list_pac_t2,Mu);

[maxef,mr,c]=EfcriticalRebarbeams(load_conditions,b,E,fdpc,arreglo_t1,arreglo_t2,...
             varDisponibles,d,h_rec,betac,disposition_rebar);

cbest=c;
maxEf=maxef;
bestMr=mr;
