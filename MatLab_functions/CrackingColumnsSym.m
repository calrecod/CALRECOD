function [InertiaXYmodif,Atransf_xy,elimxy]=CrackingColumnsSym(h,b,fdpc,rec,...
          tx,eccentricityXY,ty,Pu,cxy,conditionCrack,E)

%------------------------------------------------------------------------
% Syntax:
% [Inertia_xy_modif,Atransf_xy]=CrackingColumnsSym(h,b,fdpc,rec,tx,...
%           eccentricityXY,ty,Pu,cxy,condition_cracking,E)
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To compute the reduced inertia momentum of a column 
%          cross-section with symmetrical reinforcement considering 
%          cracking mechanisms.
% 
% OUTPUT: InertiaXYmodif: is the modified reduced inertia momentum for 
%                         both axis directions of the cross-section 
%                         considering cracking mechanisms: [Ix,Iy]
%
%         Atransf_xy:     is the transformed effective area for both axis 
%                         directions according to the cracking mechanism 
%                         (cracked or non-cracked)
%
%         elimxy:         limit eccentricity in the X and Y axis
%                         direction
%
% INPUT:  Pu:             axial load over the column's cross-section (Ton)
%
%         rec:            concrete cover for both axis direction of the 
%                         cross-section: format [coverx,covery]
%
%         fdpc:           is the factored value of f'c as 0.85f'c according 
%                         to the ACI 318-19 code  
%
%         b,h:            are the cross-section dimensions
%
%         eccentricityXY: axial load eccentricity for both axis directions:
%                         [ex,ey] 
%
%         tx,ty:          ISR width for both axis directions
%
%         cxy:            neutral axis depth for both axis directions of 
%                         the cross-section, corresponding to the optimal 
%                         reinforcement design
%
%         conditionCrack: parameter that indicates which mechanism to 
%                         consider: 
%                         format is ''Cracked'' or ''Non-cracked''
%
%         E:              Elasticity Modulus of the reinforcing steel
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

InertiaXYmodif=[];
Atransf_xy=[];
fc=fdpc/0.85;
if fc<2000 % units: kg,cm
    frot=0.8*2*sqrt(fc);
    n=E/(14000*sqrt(fc));
else       % units: lb,in
    frot=0.8*7.5*sqrt(fc);
    n=E/(57000*sqrt(fc));  
end
%% X-direction
% non-cracked cross-section inertia
rest_inertia=2*(n-1)*(b-2*rec(1))*tx*(h/2-rec(2))^2+2*(n-1)*...
                ty*(h-2*rec(2))^3; % contribution of the reinforcement
Inertia_x_modif=b*h^3/12+rest_inertia;

% transformed non-cracked section area 
a_noag=b*h+(n-1)*2*(b-2*rec(1))*tx+2*(h-2*rec(2))*(n-1)*ty;

% excentricity limit to be considered as non-cracked
elim=2*(abs(Pu)/a_noag+frot)*Inertia_x_modif/(abs(Pu)*h);
elimxy=[elim];
if conditionCrack=="Cracked"
    cracking_factor=0;
elseif conditionCrack=="Non-cracked"
    cracking_factor=inf;
end

if eccentricityXY(1)<=elim+cracking_factor
    % non-cracked cross-section
    Inertia_x_modif=Inertia_x_modif;
    Atransx=b*h+2*(b-2*rec(1))*(n-1)*tx+...
            2*(h-2*rec(2))*(n-1)*ty;
else
    % cracked cross-section:
    rest_inertia=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*tx*(h-rec(2)-cxy(1))^2+...
                (n-1)*tx*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                2*(n-1)*ty*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                2*(n-1)*ty*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                2*n*ty*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                2*n*ty*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    Inertia_x_modif=rest_inertia;
    k=cxy(1)-h/2;
    Atransx=b*(h-cxy(1))+(b-2*rec(1))*tx*((n-1)+n)+...
            2*(h/2-rec(2)-k)*ty*((n-1))+...
            2*(h/2-rec(2)+k)*ty*n;
end
InertiaXYmodif=[InertiaXYmodif,Inertia_x_modif];
Atransf_xy=[Atransf_xy,Atransx];

%% Y-direction
% non-cracked cross-section inertia
dimCol=[b,h];
b=dimCol(2); % invert cross-section dimensions
h=dimCol(1);

tvalues=[tx,ty];
tx=tvalues(2);
ty=tvalues(1);

rec_values=rec;
rec(1)=rec_values(2);
rec(2)=rec_values(1);

cvalues=cxy;
cxy=[cvalues(2),cvalues(1)];
rest_inertia=2*(n-1)*(b-2*rec(1))*tx*(h/2-rec(2))^2+2*(n-1)*...
                ty*(h-2*rec(2))^3;
Inertia_y_modif=b*h^3/12+rest_inertia;

% transformed non-cracked section area 
a_noag=b*h+(n-1)*2*(b-2*rec(1))*tx+2*(h-2*rec(2))*(n-1)*ty;

% excentricity limit to be considered as non-cracked
elim=2*(abs(Pu)/a_noag+frot)*Inertia_y_modif/(abs(Pu)*h); % cm
elimxy=[elimxy,elim];

if conditionCrack=="Cracked"
    cracking_factor=0;
elseif conditionCrack=="Non-cracked"
    cracking_factor=inf;
end

if eccentricityXY(2)<=elim+cracking_factor
    % non-cracked section:
    Inertia_y_modif=Inertia_y_modif;
    Atransy=b*(h)+2*(b-2*rec(1))*(n-1)*tx+...
            2*(h-2*rec(2))*(n-1)*ty;
else
    % cracked section:
    rest_inertia=b*cxy(1)^3/4+...
                (n)*(b-2*rec(1))*tx*(h-rec(2)-cxy(1))^2+...
                (n-1)*tx*(b-2*rec(1))*(cxy(1)-rec(2))^2+...
                2*(n-1)*ty*(h/2-rec(2)-(h/2-cxy(1)))^3/12+...
                2*(n-1)*ty*(h/2-rec(2)-(h/2-cxy(1)))*(0.5*(h/2-rec(2)-(h/2-cxy(1))))^2+...
                2*n*ty*(h/2-rec(2)+(h/2-cxy(1)))^3/12+...
                2*n*ty*(h/2-rec(2)+(h/2-cxy(1)))*(0.5*(h/2-rec(2)+(h/2-cxy(1))))^2;

    Inertia_y_modif=rest_inertia;
    k=cxy(1)-h/2;
    Atransy=b*(h-cxy(1))+(b-2*rec(1))*tx*((n-1)+n)+...
            2*(h/2-rec(2)-k)*ty*((n-1))+...
            2*(h/2-rec(2)+k)*ty*n;
end
InertiaXYmodif=[InertiaXYmodif,Inertia_y_modif];
Atransf_xy=[Atransf_xy,Atransy];
