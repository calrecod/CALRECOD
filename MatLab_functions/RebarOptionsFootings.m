function [dimb,acRebar,nv,s,arrangement]=RebarOptionsFootings(ac,dimb,...
                RebarAvailable,sepMinCode)
%------------------------------------------------------------------------
% Syntax:
% [dimb,acRebar,nv,s,arrangement]=RebarOptionsFootings(ac,dimb,...
%              RebarAvailable,sepMinCode)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal rebar option for a transversal cross
% section of a rectangular footing.
% 
% OUTPUT: dimb:                 final transversal cross-section width 
%                               dimension (in case is augmented to comply
%                               the max-min rebar area or minimum rebar 
%                               separation constraints)
%
%         acRebar:              rebar area approximate to the given ISR area
%
%         nv:                   number of rebars
%
%         arrangement:          is the rebar type to be used from the 
%                               available commercial ones
%
% INPUT:  dimb:                 transversal cross-section width
%
%         ac:                   ISR area
%
%         RebarAvailable:       commercial available rebar database: by 
%                               default size = [7,5] in format: 
%               ---------------------------------------------
%               [noption,#rebar,diam,rebarArea,lineal-weight]
%               ---------------------------------------------
%
%         sepMinCode:           minimum separation by code as 1 inch or
%                               according to the ACI 318-19 code (or 
%                               specified higher for practical construction
%                               reasons)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

ndiam=length(RebarAvailable(:,1));
matriz_best_rebar=zeros(ndiam,3);
maxiter=20;
opcion_amin=0;
j=0;
while opcion_amin==0
    j=j+1;
    
    amin_bar=inf;
    for i=1:ndiam
        trebar=i;
        dv=RebarAvailable(trebar,2); % rebar diameter
        av=(dv)^2*pi*0.25;

        s=av*dimb/ac;

        nv=fix(ac/av+1);
        if nv<2
            nv=2;
        end
        atotal=nv*av;

        matriz_best_rebar(i,1)=i;
        matriz_best_rebar(i,2)=s;
        matriz_best_rebar(i,3)=atotal;

        sepMin_02=dv;
        if sepMin_02>sepMinCode
            sepMin=sepMin_02;
        else
            sepMin=sepMinCode;
        end

        %%% discard the option that do not comply the minimum separation
        % restriction
        if s<sepMin
            matriz_best_rebar(i,:)=inf;
        end

        %%% to determine the option with the less rebar area
        if matriz_best_rebar(i,3)<amin_bar
            amin_bar=matriz_best_rebar(i,3);
            opcion_amin=i;
        end
    end
    if opcion_amin==0
        fprintf('\nDimension not fit for the footing\n');
        
        if j <= maxiter
            dimb=dimb+5.0;
        else
            acRebar=0;
            nv=0;
            s=0;
            arrangement=[];
            disp('No rebar design solution was reached for the footing!')
            break;
        end
    end
end
if opcion_amin>0
    trebar=opcion_amin;
    dv=RebarAvailable(trebar,2); % optimal rebar diameter
    av=(dv)^2*pi*0.25;

    s=av*dimb/ac; % rebar nominal separation

    nv=fix(ac/av+1); % number of rebars
    if nv<2
        nv=2;
    end
    acRebar=nv*av;

    arrangement=zeros(nv,1)+trebar;
end