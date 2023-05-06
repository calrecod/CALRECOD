function eleMec=eleMecanicosRebarColsCirc(c,fdpc,diam,rec,E,...
         rebarDisposition,av,beta1)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicosRebarColsCirc(c,fdpc,diam,rec,E,...
%        rebarDisposition,av,beta1)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the axial load and bending resistance of a 
% reinforced column of circular cross section with an ISR.
% 
% OUTPUT: eleMec:             array containing the sum of resistance forces
%                             (axial and being) as [ sum Fs,   sum Ms;
%                                                   sum Fconc, sum Mconc]
%
% INPUT:  c:                  neutral axis depth
%
%         diam:               cross-section diameter
%
%         E:                  is the Modulus of Elasticity of the 
%                             reinforcing steel
%
%         rec:                is the concrete cover
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code
%
%         beta1:              is determined as established in ACI 318 code
%                             (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-01-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

a=c/beta1;
nv=length(rebarDisposition(:,1));
eMecProfile=zeros(nv,8); % [av, xposition, yposition, d, eps, eE, Fr, Mr]

% localtion of each discrete dt over the section
for i=1:nv
    eMecProfile(i,1)=av;
    x=rebarDisposition(i,1);
    eMecProfile(i,2)=x;

    y=rebarDisposition(i,2);
    eMecProfile(i,3)=y;
end

sumaM=0;
sumaF=0;
for i=1:nv
    eMecProfile(i,4)=0.5*diam-eMecProfile(i,3); %momentum distance
    eMecProfile(i,5)=0.003/c*(eMecProfile(i,4)-c);
    
    if (eMecProfile(i,5)<-0.0021)
        eMecProfile(i,5)=-0.0021;
    elseif(eMecProfile(i,5)>0.0021)
        eMecProfile(i,5)=0.0021;
    end
    eMecProfile(i,6)=eMecProfile(i,5)*E;
    eMecProfile(i,7)=eMecProfile(i,6)*eMecProfile(i,1);
    eMecProfile(i,8)=eMecProfile(i,7)*(eMecProfile(i,4)-0.5*diam);
    
    sumaF=sumaF+eMecProfile(i,7);
    sumaM=sumaM+eMecProfile(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcretoCirc(a,fdpc,diam);

eleMec=[elemAc;
        elemConc];