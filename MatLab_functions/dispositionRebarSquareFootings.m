function [dispositionRebar]=dispositionRebarSquareFootings(b,h,rec,...
    rebarAvailable,nvt,RebarArrangement1,RebarArrangement2,axis)

%------------------------------------------------------------------------
% Syntax:
% [dispositionRebar]=dispositionRebarSquareFootings(b,h,rec,...
%   rebarAvailable,nvt,RebarArrangement1,RebarArrangement2,axis)
%
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local rebar coordinates over a transversal cross
% section of a square isolated footings for which the rebar distribution is
% placed uniformly over both transversal cross-sections.
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
%                               ______________
%                               [# rebar,diam]
%                               --------------
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
%         axis:                 is the axis direction in question: (1) axis
%                               ALONG the L footing dimension, (2) axis 
%                               ALONG the B footing dimension
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nbars=length(rebarAvailable(:,1));
if axis==1
    dextra=0;
elseif axis==2
    dextra=rec*0.7;
end
dispositionRebar=zeros(nvt(1)+nvt(2),2);
bprima=b-2*rec;
hprima=h-2*rec;

% rebar separation calculation_________________________________________

nv1=nvt(1);
nv2=nvt(2);

% rebars in tension
sumaDiam=0;
nvi=0;
for i=1:nv1
    nvi=nvi+1;
    sumaDiam=sumaDiam+rebarAvailable(RebarArrangement1(nvi),2);
end
separacion_b_inf=((bprima-sumaDiam)/(nv1-1));

% rebars in compression
sumaDiam=0;
nvi=0;
for i=1:nv2
    nvi=nvi+1;
    sumaDiam=sumaDiam+rebarAvailable(RebarArrangement2(nvi),2);
end
separation_b_sup=((bprima-sumaDiam)/(nv2-1));

%% REBAR DISTRIBUTION
% initial values for rebar disposition in tension
nivar_tension=0;
xacumInf=-bprima*0.5+0.5*rebarAvailable(nbars,2);
dispositionRebar(1,1)=xacumInf;
dispositionRebar(1,2)=-hprima*0.5+(rebarAvailable(nbars,2)*0.5)+dextra;

nivar_tension=nivar_tension+1;

% inicial values for rebar distribution in compression
nivar_comp=0;
xacumSup=-bprima*0.5+0.5*rebarAvailable(nbars,2);
dispositionRebar(nv1+1,1)=xacumSup;
dispositionRebar(nv1+1,2)=hprima*0.5-(rebarAvailable(nbars,2)*0.5)-...
                          dextra*0.5;
nivar_comp=nivar_comp+1;
   
% remaining values for coodinates........................................
% tension
for j=2:nv1
    nivar_tension=nivar_tension+1;
    xacumInf=xacumInf+0.5*rebarAvailable(RebarArrangement1(nivar_tension-1),2)+...
        separacion_b_inf+0.5*rebarAvailable(RebarArrangement1(nivar_tension),2);

    dispositionRebar(nivar_tension,1)=xacumInf;
    dispositionRebar(nivar_tension,2)=-hprima*0.5+...
                                    (rebarAvailable(nbars,2)*0.5)+dextra; 
end

% compression
for j=2:nv2
    nivar_comp=nivar_comp+1;
    xacumSup=xacumSup+0.5*rebarAvailable(RebarArrangement2(nivar_comp-1),2)+...
        separation_b_sup+0.5*rebarAvailable(RebarArrangement2(nivar_comp),2);

    dispositionRebar(nv1+nivar_comp,1)=xacumSup;
    dispositionRebar(nv1+nivar_comp,2)=hprima*0.5-...
                                rebarAvailable(nbars,2)*0.5-dextra*0.5;
end
