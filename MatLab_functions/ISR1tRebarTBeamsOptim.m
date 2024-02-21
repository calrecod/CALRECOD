function [sepbarsRestric,cbest,bestBarDisposition,bestCost,barTypesTen,...
    barTypesComp,maxEf,bestMr,areaRebar]=ISR1tRebarTBeamsOptim(bp,ht,...
    ba,ha,fc,cover,conditions,t2,pu_beams,Lb,rebarAvailable,wac)
    
%------------------------------------------------------------------------
% Syntax:
% [sepbarsRestric,cbest,b,h,bestBarDisposition,bestCost,barTypesTen,...
%  barTypesComp,maxEf,bestMr,area_var_t]=ISR1tRebarTBeamsOptim(E,bp,ht,...
%  ba,ha,fy,fc,cover,conditions,t2,pu_beams,Lb,rebarAvailable,wac)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To optimally design a rebar distribution over a beam 
% cross-section from a given reinforcement area through a T-beam the ISR.
% 
% OUTPUT: sepbarsRestric:       Is the parameter that indicates if the 
%                               minimum rebar separation restriction of 
%                               rebars in tension for the cross-section in 
%                               question is being complied: (1) indicates 
%                               that such restriction was not complied,
%                               (0) indicates that such restriction was
%                               complied
%
%         cbest:                is the depth of the neutral axis for the 
%                               optimized beam reinforced cross-section
%                               considering the optimal rebar design
%
%         bestBarDisposition:   are the local coordinates of rebar 
%                               disposition over the optimally designed 
%                               cross-section 
%
%         barTypesTen:          are the list of rebar type transformed from
%                               the ISR in tension: a vector consisting of
%                               one column of length nbars in tension
%
%         barTypesComp:         are the list of rebar type transformed from 
%                               the ISR in compression: a vector consisting 
%                               of one column of length nbars in compression
%
%         maxEf:                is the optimal final structural efficiency 
%                               for the optimal designed beam cross-section 
%                               considering the optimal rebar
%
%         bestMr:               is the bending resistance for
%                               the optimally designed T-beam cross-section 
%
%         areaRebar:            is a vector consisting of the total optimal 
%                               rebar area in tension and compression.
%                               Vector of size 1 x 2 in format:
%                               [area-tension, area-compression]
%
% INPUT:  cover:                are the concrete cover parameters 
%                               horizontally and vertically, respectively
%
%         ba:                   is the effective flange width of the T-beam 
%                               cross-section
%
%         ht:                   is total height of the T-beam cross-section
%
%         bp:                   is the web width of the T-beam 
%                               cross-section
%
%         ha:                   is the flange thickness of the T-beam
%                               cross-section
%
%         Lb:                   is the length of the beam element
%
%         fc:                   is the f'c
%
%         t2:                   is the optimal ISR consisting of vector of
%                               two elements [t_tension,t_compression]
%
%         pu_beams:             is the average unitary cost of rebar 
%                               assembly for T-beams
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

%% Determine all the possible reinforcement combinations
tma=2/3; %(in)
bpp=bp-2*cover;

ntypes=length(rebarAvailable(:,1)); % rebar diameters commercilly available
at1=t2(1)*bpp;
atmin=t2(2)*bpp;

at2=[at1,atmin];
% to determine number of rebars for each available option
areaRebar=zeros(1,2);
nvt=zeros(1,2);
for i=1:2
    ac_bar_space=zeros(ntypes,4); %[areabar separation]
    amin=inf;
    for k=1:ntypes
        j=0;
        tipovar=k;
        diamj=rebarAvailable(tipovar,2);
        areavarj=diamj^2*pi/4;
        nv=fix(at2(i)/areavarj)+1;
        if nv==1
            nv=2;
        end
        atotal=nv*areavarj;
        
        % Rebar separation:
        if fc<2000
            sepMin=max([3/2*tma*2.54,diamj]);
        else
            sepMin=max([3/2*tma,diamj]);
        end
        sep1=(bp-2*cover-nv*diamj)/(nv-1);

        sep2=(bp-2*cover-(fix(nv*0.5)+1)*diamj)/(fix(nv*0.5));

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
    areaRebar(i)=nv*best_area_var;
    nvt(i)=nv;
    if i==1
        barTypesTen=arreglo_var;
        list_pac_t1=list_pac;
    elseif i==2
        barTypesComp=arreglo_var;
        list_pac_t2=list_pac;
    end
end
Mu=conditions(1,2);

[bestBarDisposition]=RebarDisposition1tBeams(bp,ht,cover,cover,...
rebarAvailable,nvt,barTypesTen,barTypesComp,list_pac_t1,list_pac_t2,Mu);

rebarType=[barTypesTen;barTypesComp];
[maxef,mr,c]=EfRebarTBeams(conditions,bp,ht,ba,ha,Lb,fdpc,...
                           rebarType,rebarAvailable,cover,betac,...
                           bestBarDisposition);
           
cbest=c;
maxEf=maxef;
bestMr=mr;

%% Compute the construction cost
bestCost = EvaluateCostbeams(nvt(1),nvt(2),barTypesTen,barTypesComp,...
                            pu_beams,rebarAvailable,wac,Lb);
