function [disposicion_varillado]=RebarDisposition1tBeams(b,h,b_rec,h_rec,...
    varDisponibles,nv_t,arreglo_t1,arreglo_t2,list_pac_t1,list_pac_t2,Mu)

%------------------------------------------------------------------------
% Syntax:
% [disposicion_varillado]=RebarDisposition1tBeams(b,h,b_rec,h_rec,...
%   varDisponibles,nv_t,arreglo_t1,arreglo_t2,list_pac_t1,list_pac_t2,Mu)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the local coordinates of a rebar option according to
% its given data (number of rebars in tension and compression, cross-section
% dimensions, type of rebar in tension and compression and type of 
% arrangement - two pack or one pack).
% 
% OUTPUT: disposicion_varillado: local coordinates of rebars laid out over 
%                                the beam cross-section
%
% INPUT:  nv_t:                  vector containing the number of rebars in
%                                tension and compression as 
%                                [nrebar_tension,nrebar_comp]
%
%         arreglo_t1,
%         arreglo_t2:            Vectors that contain the type of rebar for
%                                the optimal option both in tension and 
%                                compression, respectively. The vectors size
%                                is of one column with nrebar rows containing
%                                a number between 1 and 7 according to the
%                                available commercial rebar types stated by
%                                default
%
%         list_pac_t1,
%         list_pac_t2:           type of arrangement for each bar, either 
%                                (1: one pack), (2: two pack). Are vectors 
%                                of one column and n-rebar rows in tension
%                                and compression, respectively
%
%         Mu:                    maximum bending load over the cross-section
%
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

disposicion_varillado=zeros(nv_t(1)+nv_t(2),2);
bprima=b-2*b_rec;
hprima=h-2*h_rec;

% To sum the rebar diameters to compute the rebar separation

% Tension
nvHor1=nv_t(1);

pacs1=length(list_pac_t1);
pacs2=length(list_pac_t2);

sumaDiam=0;
nvi=0;
for i=1:pacs1
    if list_pac_t1(i)==1
        nvi=nvi+1;
        sumaDiam=sumaDiam+varDisponibles(arreglo_t1(nvi),2);
    elseif list_pac_t1(i)==2
        nvi=nvi+2;
        sumaDiam=sumaDiam+varDisponibles(arreglo_t1(nvi-1),2);
    
    end
end
separacion_b_inf=((bprima-sumaDiam)/(pacs1-1));

% Compression
sumaDiam=0;
nvi=0;
for i=1:pacs2
    if list_pac_t2(i)==1
        nvi=nvi+1;
        sumaDiam=sumaDiam+varDisponibles(arreglo_t2(nvi),2);
    elseif list_pac_t2(i)==2
        nvi=nvi+2;
        sumaDiam=sumaDiam+varDisponibles(arreglo_t2(nvi-1),2);
    
    end
end
separacion_b_sup=((bprima-sumaDiam)/(pacs2-1));


%%--------------------------------------------------------%%%
%---------horizontal arrangement of rebars----------------%

% initial values for the rebars in tension
nivar_tension=0;
if list_pac_t1(1)==1
    xacumInf=-bprima*0.5+0.5*varDisponibles(1,2);
    disposicion_varillado(1,1)=xacumInf;
    disposicion_varillado(1,2)=-sign(Mu)*hprima*0.5+sign(Mu)*...
                                (varDisponibles(1,2)*0.5);
    
    nivar_tension=nivar_tension+1;
elseif list_pac_t1(1)==2
    xacumInf=-bprima*0.5+0.5*varDisponibles(1,2);
    disposicion_varillado(1,1)=xacumInf;
    disposicion_varillado(1,2)=-sign(Mu)*hprima*0.5+sign(Mu)*...
                                (varDisponibles(1,2)*0.5);
    
    disposicion_varillado(2,1)=xacumInf;
    disposicion_varillado(2,2)=-sign(Mu)*hprima*0.5+sign(Mu)*...
        (varDisponibles(1,2))+sign(Mu)*0.5*varDisponibles(arreglo_t1(2),2);
    
    nivar_tension=nivar_tension+2;

end

% initial values for the rebars in compression
nivar_comp=0;
if list_pac_t2(1)==1
    
    xacumSup=-bprima*0.5+0.5*varDisponibles(1,2);
    disposicion_varillado(nvHor1+1,1)=xacumSup;
    disposicion_varillado(nvHor1+1,2)=sign(Mu)*hprima*0.5-sign(Mu)*...
                                        (varDisponibles(1,2)*0.5);
    nivar_comp=nivar_comp+1;
elseif list_pac_t2(1)==2
    xacumSup=-bprima*0.5+0.5*varDisponibles(1,2);
    disposicion_varillado(nvHor1+1,1)=xacumSup;
    disposicion_varillado(nvHor1+1,2)=sign(Mu)*hprima*0.5-sign(Mu)*...
                                    (varDisponibles(1,2)*0.5);
    
    disposicion_varillado(nvHor1+2,1)=xacumSup;
    disposicion_varillado(nvHor1+2,2)=sign(Mu)*hprima*0.5-sign(Mu)*...
        (varDisponibles(1,2)*0.5)-sign(Mu)*varDisponibles(1,2);
    
    nivar_comp=nivar_comp+2;
end
   
% remaining values for coodinates........................................
for j=2:pacs1
    % coordinates of rebars in tension
    if list_pac_t1(j)==1
        nivar_tension=nivar_tension+1;
        xacumInf=xacumInf+0.5*varDisponibles(arreglo_t1(nivar_tension-1),2)+...
            separacion_b_inf+0.5*varDisponibles(arreglo_t1(nivar_tension),2);
        
        disposicion_varillado(nivar_tension,1)=xacumInf;
        disposicion_varillado(nivar_tension,2)=-sign(Mu)*hprima*0.5+...
                                    sign(Mu)*(varDisponibles(1,2)*0.5);
    elseif list_pac_t1(j)==2
        nivar_tension=nivar_tension+2;
        xacumInf=xacumInf+0.5*varDisponibles(arreglo_t1(nivar_tension-2),2)+...
            separacion_b_inf+0.5*varDisponibles(arreglo_t1(nivar_tension-1),2);
        
        disposicion_varillado(nivar_tension-1,1)=xacumInf;
        disposicion_varillado(nivar_tension-1,2)=sign(Mu)*(-hprima*0.5+...
                                            (varDisponibles(1,2)*0.5));
        
        disposicion_varillado(nivar_tension,1)=xacumInf;
        disposicion_varillado(nivar_tension,2)=sign(Mu)*(-hprima*0.5+...
                        (varDisponibles(1,2))+0.5*varDisponibles(1,2));
        
    end    
end

% rebar in compression
for j=2:pacs2
    if list_pac_t2(j)==1
        nivar_comp=nivar_comp+1;
        % coordinates of rebars in compression
        xacumSup=xacumSup+0.5*varDisponibles(arreglo_t2(nivar_comp-1),2)+separacion_b_sup+...
            0.5*varDisponibles(arreglo_t2(nivar_comp),2);
        
        disposicion_varillado(nvHor1+nivar_comp,1)=xacumSup;
        disposicion_varillado(nvHor1+nivar_comp,2)=sign(Mu)*(hprima*0.5-...
                                                   varDisponibles(1,2)*0.5);
        
    elseif list_pac_t2(j)==2
        
        nivar_comp=nivar_comp+2;
        
        xacumSup=xacumSup+0.5*varDisponibles(arreglo_t2(nivar_comp-1),2)+...
            separacion_b_sup+0.5*varDisponibles(arreglo_t2(nivar_comp-1),2);
        
        disposicion_varillado(nvHor1+nivar_comp-1,1)=xacumSup;
        disposicion_varillado(nvHor1+nivar_comp-1,2)=sign(Mu)*(hprima*0.5-...
                                                varDisponibles(1,2)*0.5);
        
        disposicion_varillado(nvHor1+nivar_comp,1)=xacumSup;
        disposicion_varillado(nvHor1+nivar_comp,2)=sign(Mu)*(hprima*0.5-...
                            varDisponibles(1,2)-varDisponibles(1,2)*0.5);
    end
end
