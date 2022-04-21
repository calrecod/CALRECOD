function [d,sigma13,sigma24,sigma12,sigma34,a13,a24,a12,a34,qmax13,...
    qmax24,qmax12,qmax34]=shearFootings(be,le,qprom,dimCol,qu01,qu02,qu03,...
    qu04,pu,d,fc)

%------------------------------------------------------------------------
% Syntax:
% [d,sigma13,sigma24,sigma12,sigma34,a13,a24,a12,a34,qmax13,...
%    qmax24,qmax12,qmax34]=shearFootings(be,le,qprom,dimCol,qu01,qu02,qu03,...
%    qu04,pu,d,fc)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the demand of shear stresses and resistant ones of an
% isolated footing subject to biaxial eccentric actions, considering two 
% mechanisms (punching and beam shearing).
% 
% OUTPUT: d:               effective modified height dimension based on the
%                          critical acting shear demand over the isolated 
%                          footing
%
%         a13,a24,a12,a34: are max effective unit-stress acting on each axis 
%                          direction of the footing due to the contact 
%                          pressures (see documentaiton) for each of its 
%                          four plan boundaries
%
%         qmax13,qmax24
%         qmax12,qmax34:   are the max contact pressures for each of the 
%                          four footing boundaries at their respective ends
%                          (see documentation)
%
%         sigma13,sigma24
%         sigma12,sigma34: define the stress magnitude at the other end of
%                          the block unit stress (see documentaiton) for 
%                          each of the four footing boundaries
%
% INPUT:  qprom:           is the average contact pressure considering the
%                          pressures at the four footing's corners (Kg/cm2)
%
%         dimCol:          are the cross-section dimensions of the column
%                          that the footing supports (cm)
%
%         qu01,qu02
%         qu03,qu04:       are the contact pressure magnitudes at the four 
%                          footing's corners (see documentation)
%
%         be,le:           are the transversal dimensions of the isolated 
%                          footing on plan view (cm)
%
%         pu:              is the axial load reaction from the column (Kg)
%
%         d:               is the effective footing height (cm)
%
%         fc:              is the f'c used for the footing (Kg/cm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%%%%%%%%% Shear revision _______________________________________________

%%%%%%%% Punching shear %%%%%%%%%
vneta=pu-qprom*dimCol(1)*dimCol(2);
%vneta=q_prom*(be*le-((dimCol(1)+d)*(dimCol(2)+d)))
acorte=d*((dimCol(1)+d)*2+(dimCol(2)+d)*2);

vup=vneta/acorte; % punching shear demand
vcp=0.85*fc^0.5; % resistant punching shear

dpunz=d;
if vup>vcp
    while vup>vcp
        dpunz=dpunz+5; % minimum required
        acorte=dpunz*((dimCol(1)+dpunz)*2+(dimCol(2)+dpunz)*2);

        vup=vneta/acorte;

    end
elseif vup<vcp
    while vup<vcp
        dpunz=dpunz-5; % minimum required
        acorte=dpunz*((dimCol(1)+dpunz)*2+(dimCol(2)+dpunz)*2);

        vup=vneta/acorte;
        if dpunz<15
            break;
        end
    end
    dpunz=dpunz+5; % to adjust to the limit (minimum height required ...
                   % punching shear demand
end
%%%%%%%%%%%% Revision as a beam ________________________________________

%%%% it is considered that (1-3) y (2-4) go along the Le dimension %%%%%

if qu01>qu03
    qmax13=qu01;
    qmin13=qu03;
else
    qmax13=qu03;
    qmin13=qu01;
end

ma13=(qmax13-qmin13)/be;

if qu02>qu04
    qmax24=qu02;
    qmin_24=qu04;
else
    qmax24=qu04;
    qmin_24=qu02;
end

ma24=(qmax24-qmin_24)/be;

%%%% se considera que 3-4 y 1-2 van a lo largo de Be %%%%%

if qu03>qu04
    qmax34=qu03;
    qmin34=qu04;
else
    qmax34=qu04;
    qmin34=qu03;
end

ma34=(qmax34-qmin34)/le;

if qu01>qu02
    qmax12=qu01;
    qmin12=qu02;
else
    qmax12=qu02;
    qmin12=qu01;
end

ma12=(qmax12-qmin12)/le;

%%%%%%%%%%%%%% sentido de be %%%%%%%%%%%%%%%%

bp=(le-dimCol(1)-d)*0.5;
sigma13=qmax13-ma13*bp;

a13=sigma13*bp+(qmax13-sigma13)*bp*1/2;

sigma24=qmax24-ma24*bp;
a24=sigma24*bp+(qmax24-sigma24)*bp*1/2;

aprom=(a13+a24)*0.5;
fvi=aprom*be;

vv_be=fvi/(be*d); % demand of shear as a beam
vcv_be=0.5*0.85*fc^0.5; % resistant shear as a beam

dpbe=fvi/(be*vcv_be); % minimum height required as a beam

%%%%%%%%%%%%%% on the le-direction %%%%%%%%%%%%%%%%

lp=(be-dimCol(2)-d)*0.5;
sigma12=qmax12-ma12*lp;
a12=(sigma12*lp+((qmax12-sigma12)*lp*1/2));

sigma34=qmax34-ma34*lp;
a34=(sigma34*lp+((qmax34-sigma34)*lp*0.5));

aprom=(a12+a34)*0.5;
fvi=aprom*le;

vv_le=fvi/(le*d); % demand of shear as a beam
vcv_le=0.5*0.85*fc^0.5; % resistant shear as a beam

dple=fvi/(le*vcv_le); % minimum height required as a beam

%%% to determined the max required height considering both mechanisms
vector_d=[dpunz,dpbe,dple];
d=max(vector_d);

d=d+5-mod(d,5); % effectve modified height
