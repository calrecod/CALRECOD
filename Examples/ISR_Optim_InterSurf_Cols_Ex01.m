% ISR_Optim_InterSurf_Cols_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    1.- To design determine the optimal reinforcement area of a
%        rectangular column cross-section subject to biaxial bending
%        through the computation of interaction surfaces.
%
%    2.- To optimally design the rebar on the column cross-section 
%        once an approximate optimal reinforcement area has been found.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%----------------------------------------------------------------

clear all
clc

%% GEOMETRY
height=300; % column's length
b=30; % cross-section width
h=30; % cross-section height

%% MATERIALS
fy=4200; % yield stress of rebars 
fc=250; % concrete's compressive strength Kg/cm2
ws=0.0078; % volumetric weight of reinforcing steel (Kg/cm3)

if fc<280
    beta1=0.85;
elseif fc>=280
    beta1=1.05-fc/1400;
    if (beta1<0.65)
        beta1=0.65;
    elseif (beta1>0.85)
        beta1=0.85;
    end
end

%% Loads
load_conditions=[1 -53.52e3 -9.602e5 10e5;
                 2 -40.22e3 16.4e5 -9e5]; % [nload, Pu, Mx, My]

%% ADITIONAL PARAMETERS
rec=[4 4]; % [coverx covery] - concrete cover

cols_sym_asym_isr="Symmetric";
if cols_sym_asym_isr=="Symmetric" || cols_sym_asym_isr=="Asymmetric"
    pu_cols=[29.19, 29.06, 28.93, 28.93, 28.93, 28.93, 28.93];

elseif  cols_sym_asym_isr=="ISR"
    pu_cols=[28.93];
end

condition_cracking="Cracked"; % "Cracked" or "Non-cracked"

%% Rebar data
% Commerically available rebars
RebarAvailable=[4 4/8*2.54;
                5 5/8*2.54;
                6 6/8*2.54;
                8 8/8*2.54;
                9 9/8*2.54;
                10 10/8*2.54;
                12 12/8*2.54];
            
% Ductility demand
ductility=3;

%% ISR Optimization
[b,h,cost_elem_col,AsISR,Ef_sec_col,Mr_col,t_value_x,...
t_value_y,cxp,bestLoad0]=isrColsInterSurf(pu_cols,height,b,h,rec,fy,fc,...
load_conditions,0.0078,ductility,1,1);

%% Optimization rebar design
puCostCardBuild=[1.2,0.9,128,214.75,0.04,0.105,7];

Es=fy/0.0021; % Moudulus of Elasticity of reinforcing steel
npdiag=30; % number of points to be computed for the interaction diagrams
OptionRP=10; % select the rebar prototype for optimal design
rangeCFA=[0,1]; % Acceptable range of constructability for the rebar designs
if OptionRP==1
    %% Asymmetrical rebar design
    % RP: (Asym-2pack-Sym4Diam)

    [Mr_col,h,bestArea1,bestCost,bestdiagram,...
    bestnv,bestEf,bestArrangement,bestDisposition,nv4,bestcxy,bestCP,...
    bestLoad,bestCFA1]=Asym2packSym4DiamIntSurf(b,h,rec,AsISR,Es,npdiag,fc*0.85,...
    beta1,load_conditions,ws,height,RebarAvailable,ductility,...
    puCostCardBuild,rangeCFA);

elseif OptionRP==2
    % RP: (Asym-2pack-1Diam)
    
    [bestMr,h,bestArea2,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestCP,bestLoad2,bestCFA2]=Asym2pack1DiamIntSurf(b,...
    h,rec,AsISR,npdiag,fc*0.85,fy,beta1,load_conditions,RebarAvailable,...
    height,ws,pu_asym_cols,ductility,puCostCardBuild,rangeCFA);

elseif OptionRP==3
    % RP: (Asym-2pack-4Diam)
    
    [Mr_col,h,bestArea3,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,estcxy,bestCP,bestLoad3,bestCFA3]=Asym2pack4DiamIntSurf(b,h,...
    rec,AsISR,npdiag,fc*0.85,beta1,fy,load_conditions,pu_asym_cols,...
    RebarAvailable,ws,height,ductility,puCostCardBuild,rangeCFA);

elseif OptionRP==4
    % RP: (Asym-Sym-4Diam)
    
    [Mr_col,h,bestArea4,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestCP,bestLoad4,bestCFA4]=AsymmSym4DiamIntSurf(b,h,...
    rec,AsISR,Es,npdiag,fc*0.85,beta1,height,ws,RebarAvailable,...
    load_conditions,pu_asym_cols,ductility,puCostCardBuild,rangeCFA);

elseif OptionRP==5
    % RP: (Asym-4Diam)
    
    [Mr_col,h,bestArea5,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestCP,bestLoad5,bestCFA5]=Asym4DiamIntSurf(b,h,rec,...
    AsISR,Es,npdiag,fc*0.85,beta1,RebarAvailable,height,ws,load_conditions,...
    pu_asym_cols,ductility,puCostCardBuild,rangeCFA);

elseif OptionRP==6
    % RP: (Asym-1Diam)
    pu_asym_cols=1/0.75*sum(pu_cols)/length(pu_cols); % asymmetrical
                                                     % rebar
    [Mr_col,h,bestArea6,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestCP,bestLoad6,bestCFA6]=Asym1DiamIntSurf(b,h,rec,...
    AsISR,npdiag,fc*0.85,beta1,load_conditions,pu_asym_cols,height,ws,...
    RebarAvailable,ductility,puCostCardBuild,rangeCFA);

elseif OptionRP<=6 && OptionRP>=1
    section=[0.5*b 0.5*h;
             -0.5*b 0.5*h;
             -0.5*b -0.5*h;
              0.5*b -0.5*h;
              0.5*b 0.5*h];

    PlotRotRecSecRebarCols(section,bestLoad6,bestdiagram,...
                        bestDisposition,b,h,bestArrangement);
                    
elseif OptionRP==7
    %% Symmetrical rebar design
    % RP: (Sym-12Diam)
    
    [Mr_col,h,bestArea7,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestLoad7,bestCFA7]=supOptimRebarSymIntSurf...
    (b,h,rec,AsISR,Es,npdiag,fc*0.85,beta1,load_conditions,ws,height,...
    ductility,RebarAvailable,puCostCardBuild,rangeCFA,1);

elseif OptionRP==8
    % RP: (Sym-1Diam)
    
    [Mr_col,h,bestArea8,lowestCost,ovMostEc,nvEc,maxEfEc,bestArrangement,...
    bestDisposition,nvxy,bestdiagram,bestLoad8,bestCFA8]=optrebarColsSymIntSurf(b,h,...
    rec,AsISR,Es,npdiag,fc*0.85,beta1,RebarAvailable,ws,height,...
    load_conditions,puCostCardBuild,rangeCFA,1);

elseif OptionRP==9
    % RP: (Sym-2pack-1Diam)
    
    [Mr_col,h,bestArea9,lowestCost,ovMostEc,nvEc,maxEfEc,bestArrangement,...
    bestDisposition,nvxy,bestdiagram,bestLoad9,bestCFA9]=optRebarColsSym2PackIntSurf...
    (b,h,rec,AsISR,Es,npdiag,fc*0.85,beta1,RebarAvailable,ws,...
    height,load_conditions,puCostCardBuild,rangeCFA,1);

elseif OptionRP==10
    % RP: (Sym-2pack-2Diam)
    
    [Mr_col,h,bestArea10,bestCost,bestdiagram,bestnv,bestEf,bestArrangement,...
    bestDisposition,nv4,bestcxy,bestLoad10,bestCFA10]=supOptRebarSym2PackIntSurf(b,h,...
    rec,AsISR,Es,npdiag,fc*0.85,beta1,load_conditions,ws,height,...
    RebarAvailable,ductility,puCostCardBuild,rangeCFA,1);

end
