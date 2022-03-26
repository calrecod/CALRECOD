function [InertiaXYmodif,Atransfxy]=CrackingColumnsAsym(h,b,fdpc,rec,eccentricityXY,...
                            t1bar,t2bar,t3bar,t4bar,Pu,cxy,conditionCracking,cp)

%------------------------------------------------------------------------
% Syntax:
% [InertiaXYmodif,Atransfxy]=CrackingColumnsAsym(h,b,fdpc,rec,eccentricityXY,...
%                            t1bar,t2bar,t3bar,t4bar,Pu,cxy,conditionCracking,cp)
%------------------------------------------------------------------------
% PURPOSE: To compute the reduced inertia momentum of a column 
%          cross-section with asymmetrical reinforcement considering 
%          cracking mechanisms.
% 
% OUTPUT: InertiaXYmodif:   is the modified reduced inertia momentum for 
%                           both axis directions of the cross-section 
%                           considering cracking mechanisms: [Ix,Iy]
%
%         Atransfxy:        is the transformed effective area for both axis 
%                           directions according to the cracking mechanism 
%                           (cracked or non-cracked)
%
% INPUT:  Pu:               axial load over the column's cross-section
%
%         rec:              concrete cover for both axis direction of the 
%                           cross-section: format [coverx,covery]
%
%         fdpc:             is the factored value of f'c as 0.85f'c according 
%                           to the ACI 318-19 code  
%
%         b,h:              are the cross-section dimensions
%
%         eccentricityXY:   axial load eccentricity for both axis directions:
%                           [ex,ey]
%
%         t1bar,t2bar:      ISR horizontal width (X-axis)
%
%         t3bar,t4bar:      ISR vertical width (Y-axis)
%
%         cxy:              neutral axis depth for both axis directions of 
%                           the cross-section, corresponding to the optimal 
%                           reinforcement design
%
%         conditionCrack:   parameter that indicates which mechanism to 
%                           consider neglecting the rupture modulus fr_ot: 
%                           format is ''Cracked'' or ''Non-cracked''
%
%         cp:               Plastic Center depth (cm)
%
%         E:                Elasticity Modulus (Kg/cm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

frot=0.8*2*sqrt(fdpc);
InertiaXYmodif=[];
Atransfxy=[];
E=2.0e6;

n=E/(14000*sqrt(fdpc));

% Inertia X-axis_____________________________________________________
t1_var_der=t2bar;
t2_var_der=t1bar;
t3_var_der=t4bar;
t4_var_der=t3bar;

% non-cracked cross-section inertia
rest_inertia=b*h*((h/2-cp(1))^2)+...
            (n-1)*(b-2*rec(1))*(t1bar)*((rec(2)-cp(1))^2)+...
            (n-1)*t2bar*(b-2*rec(1))*((h-rec(2)-cp(1))^2)+...
            (n-1)/12*(t3bar+t4bar)*((h-2*rec(2))^3)+...
            (n-1)*(t3bar+t4bar)*(h-2*rec(2))*((h/2-cp(1))^2);

Inertia_x_modif=b*h^3/12+rest_inertia;

% transformed non-cracked section area 
a_noag=b*h+(n-1)*(b-2*rec(1))*(t1bar+t2bar)+(h-2*rec(2))*(n-1)*(t3bar+t4bar);

% excentricity limit to be considered as non-cracked
elim=2*(abs(Pu*1000)/a_noag+frot)*Inertia_x_modif/(abs(Pu*1000)*h);

if conditionCracking=="Cracked"
    cracking_factor=0;
elseif conditionCracking=="Non-cracked"
    cracking_factor=inf;
end

if eccentricityXY(1)<=elim+cracking_factor
    % sección no-agrietada
    Inertia_x_modif=Inertia_x_modif;
    Atransx=b*h+(b-2*rec(1))*(n-1)*t1bar+...
            (b-2*rec(1))*(n-1)*t2bar+...
            (h-2*rec(2))*(n-1)*t3bar+...
            (h-2*rec(2))*(n-1)*t4bar;
else
    
    % cracked section ---- (compression on the upper part of the section)
    rest_inertia=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*t2bar*(h-rec(2)-cxy(1))^2+...
                (n-1)*t1bar*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    % cracked section ---- (compression on the lower part of the section)
    cxy_right=cxy; % The same neutral axis depth value is considered for...
                   % both directions of analysis.
    rest_inertia_der=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*t1_var_der*(h-rec(2)-cxy(1))^2+...
                (n-1)*t2_var_der*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    % THE LOWER INERTIA IS USED (FOR BEING THE CRITICAL ONE)...
    if rest_inertia_der<rest_inertia
        rest_inertia=rest_inertia_der;
    end
    k=cxy(1)-h/2;
    Inertia_x_modif=b*h^3/12+rest_inertia;
    Atransx=b*(h-cxy(1))+(b-2*rec(1))*t1bar*((n-1))+...
            +(b-2*rec(1))*t2bar*((n))+...
            2*(h/2-rec(2)-k)*t3bar*((n-1))+...
            2*(h/2-rec(2)+k)*t4bar*n;
end
InertiaXYmodif=[InertiaXYmodif,Inertia_x_modif];
Atransfxy=[Atransfxy,Atransx];

% Y-direction_________________________________________________________

cp_values=cp;
cp(1)=cp_values(2);
cp(2)=cp_values(1);
dimCol=[b,h];
b=dimCol(2);
h=dimCol(1);

tvalues=[t1bar,t2bar,t3bar,t4bar];
t1bar=tvalues(4);
t2bar=tvalues(3);
t3bar=tvalues(1);
t4bar=tvalues(2);

t1_var_der=t2bar;
t2_var_der=t1bar;
t3_var_der=t4bar;
t4_var_der=t3bar;

rec_values=rec;
rec(1)=rec_values(2);
rec(2)=rec_values(1);

cvalues=cxy;
cxy=[cvalues(2),cvalues(1)];

% non-cracked cross-section inertia
rest_inertia=b*h*((h/2-cp(1))^2)+...
            (n-1)*(b-2*rec(1))*(t1bar)*((rec(2)-cp(1))^2)+...
            (n-1)*t2bar*(b-2*rec(1))*((h-rec(2)-cp(1))^2)+...
            (n-1)/12*(t3bar+t4bar)*((h-2*rec(2))^3)+...
            (n-1)*(t3bar+t4bar)*(h-2*rec(2))*((h/2-cp(1))^2);

Inertia_y_modif=b*h^3/12+rest_inertia;

% transformed non-cracked section area 
a_noag=b*h+(b-2*rec(1))*(n-1)*t1bar+...
            b*h+(b-2*rec(1))*(n-1)*t2bar+...
            (h-2*rec(2))*(n-1)*t3bar+...
            (h-2*rec(2))*(n-1)*t4bar;
% excentricity limit to be considered as non-cracked
elim=2*(abs(Pu*1000)/a_noag+frot)*Inertia_y_modif/(abs(Pu*1000)*h);

if conditionCracking=="Cracked"
    cracking_factor=0;
elseif conditionCracking=="Non-cracked"
    cracking_factor=inf;
end

if eccentricityXY(2)<=elim+cracking_factor
    % Non-cracked section
    Inertia_y_modif=Inertia_y_modif;
    Atransy=b*h+(b-2*rec(1))*(n-1)*t1bar+...
            (b-2*rec(1))*(n-1)*t2bar+...
            (h-2*rec(2))*(n-1)*t3bar+...
            (h-2*rec(2))*(n-1)*t4bar;
else

    % cracked section ---- (compression on the upper part of the section)
    rest_inertia=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*t2bar*(h-rec(2)-cxy(1))^2+...
                (n-1)*t1bar*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    % cracked section ---- (compression on the lower part of the section)
    cxy_right=cxy; % The same neutral axis depth value is considered for...
                   % both directions of analysis.
    rest_inertia_der=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*t1_var_der*(h-rec(2)-cxy(1))^2+...
                (n-1)*t2_var_der*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                (n-1)*(t3bar+t4bar)*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                n*(t3bar+t4bar)*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    %  THE LOWER INERTIA IS USED (FOR BEING THE CRITICAL ONE)...
    if rest_inertia_der<rest_inertia
        rest_inertia=rest_inertia_der;
    end
    Inertia_y_modif=b*h^3/12+rest_inertia;
    k=cxy(1)-h/2;
    Atransy=b*(h-cxy(1))+(b-2*rec(1))*t1bar*((n-1))+...
            (b-2*rec(1))*t2bar*n+...
            2*(h/2-rec(2)-k)*t3bar*((n-1))+...
            2*(h/2-rec(2)+k)*t4bar*n;
end
InertiaXYmodif=[InertiaXYmodif,Inertia_y_modif];
Atransfxy=[Atransfxy,Atransy];

