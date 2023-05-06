function [dispositionRebar,arrangementfinal]=dispositionRebarRectangularFootings(b,h,rec,...
    rebarAvailable,nvt,arrangement1,arrangement2,axis,largerDim,dim_zap)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar,arrangementfinal]=dispositionRebarRectangularFootings(b,h,rec,...
%   rebarAvailable,nvt,arrangement1,arrangement2,axis,largerDim,dim_zap)
%
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local rebar coordinates over a transversal cross
% section of a rectangular isolated footing for which the rebar distribution
% is placed non-uniformly over the larger transversal cross-sections.
% 
% OUTPUT: dispositionRebar:     is the array containing the local rebar 
%                               coordinates over the transversal cross
%                               section
%
% INPUT:  b,h:                  are the transversal cross-section dimensions
%
%         rec:                  is the concrete cover
%
%         rebarAvailable:       are the commercial available rebar database:
%                               size = n# x 2, by default, in format: 
%                           ______________
%                           [# rebar,diam]
%                           --------------
%
%         nvt:                  is a vector containing the quantity of rebars
%                               both in compression and tension as 
%                   _______________________________
%                   [nb_{tension},nb_{compression}]
%                   -------------------------------
%
%         RebarArrangement1,
%         RebarArrangement2:    are the arrays containing the rebar type 
%                               used both in tension and compression, 
%                               respectively. Size = [nbars,1] in format: 
%                   __________________________
%                   [#rebar{1},...,#rebar_{n}]
%                   --------------------------
%
%         axis:                 is the axis direction in question: (1) axis
%                               ALONG the L footing dimension, (2) axis 
%                               ALONG the B footing dimension
%
%         largerDim:            is the parameter that indicates which 
%                               transversal cross-section is in analysis: 
%                               for rectangular footings (1) means that the
%                               larger dimension is the L dimensions, (2)
%                               means that the larger dimension is the B
%                               dimension. If the footing is squared, then
%                               the parameters takes the value of (3)
%
%         dim_zap:              are the two transversal cross-section
%                               dimensions of the footing: in format [B,L]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
ndiam=length(rebarAvailable(:,1));
if axis==1
    dextra=0;
elseif axis==2
    dextra=rec*0.7;
end

disposicion_varillado=zeros(nvt(1)+nvt(2),2);
    
bz=dim_zap(1);
lz=dim_zap(2);

% To determine which dimensions width is larger and which is shorter
if largerDim==1
    acprima=nvt(1)*rebarAvailable(arrangement1(1),2)^2*pi/4*(2*bz/(lz+bz));
    lextremos=(lz-bz)*0.5;
    bzp=bz;
elseif largerDim==2
    
    acprima=nvt(1)*rebarAvailable(arrangement1(1),2)^2*pi/4*(2*lz/(lz+bz));
    lextremos=(bz-lz)*0.5;
    bzp=lz;
end

%%%%%%%%% Determine number of rebars for the reduced dimension %%%%%%%%%
av_red=0;
nv_red=0;
while av_red<acprima
    nv_red=nv_red+1;
    av_red=av_red+rebarAvailable(arrangement1(1),2);
end

nv_nonred=length(arrangement1);
if mod(nv_nonred,2)~=0 % if the number of rebars is impar
    if nv_red<nv_nonred
        if mod(nv_red,2)==0
            nv_red=nv_red+1;
        end
    end
elseif mod(nv_nonred,2)==0 % if the number of rebars is par
    if nv_red<nv_nonred
        if mod(nv_red,2)~=0
            nv_red=nv_red+1;
        end
    end
end
ac_rest=(nvt(1)-nv_red)*rebarAvailable(arrangement1(1),2);

nv_ends=(nvt(1)-nv_red)*0.5;

%%%%% To determine if it is necessary to distribute rebars in packs 
%%%%% at the ends of the larger footing dimension______________________
nv_pacs=0;
list_pacs1_red=[];
rem_bars=nv_red-nv_pacs;
while nv_pacs<nv_red
    if rem_bars>1
        nv_pacs=nv_pacs+1;
        rem_bars=nv_red-nv_pacs;
        list_pacs1_red=[list_pacs1_red,1];
    elseif rem_bars==1
        nv_pacs=nv_pacs+1;
        rem_bars=nv_red-nv_pacs;
        list_pacs1_red=[list_pacs1_red,1];
    end
end
hprima=h-2*rec;
if nv_pacs==nv_red
    list_pacs1_red=zeros(nv_pacs)+1;
end

% rebar separation calculation_____________________________________

% number of rebars in tension

if nv_red<nvt(1)
    nv1=nv_red;
    bprima=bzp;
elseif nv_red>=nvt(1)
    nv1=nvt(1);
    bprima=b-2*rec;
end
nv2=nvt(2); % number of rebars in compression 

%%%%%%%%%%%%% rebar separation calculation in tension for the 
% reduced dimension _________________________________________________
sumaDiam=0;
nvi=0;
for i=1:nv1
    nvi=nvi+1;
    sumaDiam=sumaDiam+rebarAvailable(arrangement1(nvi),2);
end
separacion_b_inf=((bprima-sumaDiam)/(nv1-1));

% REBAR DISTRIBUTION (Uniformly in individual rebars)__________________
% initial values of rebar coordinates in tension
nivar_tension=0;
xacumInf=-bprima*0.5+0.5*rebarAvailable(ndiam,2);
disposicion_varillado(1,1)=xacumInf;
disposicion_varillado(1,2)=-hprima*0.5+(rebarAvailable(ndiam,2)*0.5)+dextra;

nivar_tension=nivar_tension+1;

% remaining values for coodinates........................................
% tension
for j=2:nv1
    nivar_tension=nivar_tension+1;
    xacumInf=xacumInf+0.5*rebarAvailable(arrangement1(nivar_tension-1),2)+...
        separacion_b_inf+0.5*rebarAvailable(arrangement1(nivar_tension),2);

    disposicion_varillado(nivar_tension,1)=xacumInf;
    disposicion_varillado(nivar_tension,2)=-hprima*0.5+...
                                    (rebarAvailable(ndiam,2)*0.5)+dextra;
end

bprima=b-2*rec;
% compresion
sumaDiam=0;
nvi=0;
for i=1:nv2
    nvi=nvi+1;
    sumaDiam=sumaDiam+rebarAvailable(arrangement2(nvi),2);
end
separacion_b_sup=((bprima-sumaDiam)/(nv2-1));

% initial values of rebar coordinates in compression
nivar_comp=0;
xacumSup=-bprima*0.5+0.5*rebarAvailable(ndiam,2);
disposicion_varillado(nv1+1,1)=xacumSup;
disposicion_varillado(nv1+1,2)=hprima*0.5-(rebarAvailable(ndiam,2)*0.5)-...
                            dextra*0.5;
nivar_comp=nivar_comp+1;

% remaining values for coodinates........................................
% compression
for j=2:nv2
    nivar_comp=nivar_comp+1;
    xacumSup=xacumSup+0.5*rebarAvailable(arrangement2(nivar_comp-1),2)+...
        separacion_b_sup+0.5*rebarAvailable(arrangement2(nivar_comp),2);

    disposicion_varillado(nv1+nivar_comp,1)=xacumSup;
    disposicion_varillado(nv1+nivar_comp,2)=hprima*0.5-...
                                rebarAvailable(ndiam,2)*0.5-dextra*0.5;

end

%%%%%%% REBAR DISTRIBUTION AT THE ENDS OF THE LARGER FOOTING DIMENSION____
%%%%%%% rebar separation calculation %%%%%%%%%%%%%%

%%%%% crear otra lista de pacs en tension %%%%%
sumaDiam=nv_ends*rebarAvailable(arrangement1(1),2);
tma=0.75; %pulg
sepMin=3/2*tma*2.54;

nv_pacs_ext=0;
list_pacs1_ext=[];
rem_bars_ext=nv_ends-nv_pacs_ext;
while nv_pacs_ext<nv_ends
    nv_pacs_ext=nv_pacs_ext+1;
    list_pacs1_ext=[list_pacs1_ext,1];
end
npacs_ext=length(list_pacs1_ext);
separacion_ext=(((lextremos-rec-separacion_b_inf)-sumaDiam)/(npacs_ext-1));

bprima=b-2*rec;
% REBAR DISTRIBUTION ___________________________________________________
if nv_ends~=0
    % initial coordinate values for rebars in tension ( left end )
    nivar_tension=0;
    xacumInf=-bprima*0.5+0.5*rebarAvailable(ndiam,2);
    disposicion_varillado(nv_red+1,1)=xacumInf;
    disposicion_varillado(nv_red+1,2)=-hprima*0.5+...
                                    (rebarAvailable(ndiam,2)*0.5)+dextra;

    nivar_tension=nivar_tension+1;

    % remaining values for coodinates ( left end ) ...............
    for j=2:npacs_ext
        
        nivar_tension=nivar_tension+1;
        xacumInf=xacumInf+0.5*rebarAvailable(arrangement1(nivar_tension-1),2)+...
            separacion_ext+0.5*rebarAvailable(arrangement1(nivar_tension),2);

        disposicion_varillado(nv_red+nivar_tension,1)=xacumInf;
        disposicion_varillado(nv_red+nivar_tension,2)=-hprima*0.5+...
                                        (rebarAvailable(ndiam,2)*0.5)+dextra;
   
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% REBAR DISTRIBUTION AT THE RIGHT END ...............-----------

    % Initial coordinate values of rebars in tension ( right end )

    xacumInf=bprima*0.5-0.5*rebarAvailable(ndiam,2);
    disposicion_varillado(nv_red+nivar_tension+1,1)=xacumInf;
    disposicion_varillado(nv_red+nivar_tension+1,2)=-hprima*0.5+...
                                    (rebarAvailable(ndiam,2)*0.5)+dextra;

    nivar_tension=nivar_tension+1;
    
    % remaining values for coodinates of rebars in tension ( right end )
    for j=2:npacs_ext
        nivar_tension=nivar_tension+1;
        xacumInf=xacumInf-0.5*rebarAvailable(arrangement1(nivar_tension-1),2)-...
            separacion_ext-0.5*rebarAvailable(arrangement1(nivar_tension),2);

        disposicion_varillado(nv_red+nivar_tension,1)=xacumInf;
        disposicion_varillado(nv_red+nivar_tension,2)=-hprima*0.5+...
                                      (rebarAvailable(ndiam,2)*0.5)+dextra;
  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%-----------------------%%%%%%%%%%%%%%%%%%%%%%%%

dispositionRebar=disposicion_varillado;
arrangementfinal=arrangement1;
