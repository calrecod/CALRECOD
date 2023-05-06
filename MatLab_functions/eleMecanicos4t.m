function eleMec=eleMecanicos4t(c,a,fdpc,h,b,rec,E,t1,t2,t3,t4)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicos4t(c,a,fdpc,h,b,rec,E,t1,t2,t3,t4)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To determine the axial load and bending resistance of a 
% rectangular reinforced column cross section with an ISR for a correspon-
% ding neutral axis value c. A discretization approach is used for the ISR.
% 
% OUTPUT: eleMec:             array containing the sum of resistance forces
%                             (axial and being) as [ sum Fs,   sum Ms;
%                                                   sum Fconc, sum Mconc]
%
% INPUT:  b,h:                cross-section dimensions of column (width 
%                             and height)
%
%         t1,t2,t3,t4:        are the ISR widths of a 4t-ISR (strictly for
%                             columns)
%
%         rec:                is a vector containing the concrete cover for
%                             both cross-section axis as [coverX,coverY]
%
%         fdpc:               is the reduced value of f'c with the factor 
%                             0.85 as prescribed in the ACI 318-19 code:
%                             fdpc=0.85*fc
%                               
%         beta1:              is determined as established in ACI 318 code
%                             (see Documentation)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bp=b-2*rec(1);
hp=h-2*rec(2);

ndt=100;
eMecProfile=zeros(ndt*4,8);

% [dA, xposition, yposition, d, eps, eE, Fr, Mr]

%calculate dA for each t area
dA1=(bp)*t1/ndt;
dA2=(bp)*t2/ndt;
dA3=(hp)*t3/ndt;
dA4=(hp)*t4/ndt;

% localtion of each discrete dt over the section
for i=1:ndt*4
    if i<=ndt*1
        j=i; %counting the number of dt along a face
        eMecProfile(i,1)=dA1;
        x=-0.5*b+rec(1)+(j-1)*bp/ndt+0.5*(bp/ndt);
        eMecProfile(i,2)=x;
        
        y=0.5*h-rec(2)-0.5*t1;
        eMecProfile(i,3)=y;
        
    elseif i>ndt*1 && i<=ndt*2
        j=i-ndt; %counting the number of dt along a face
        eMecProfile(i,1)=dA2;
        x=-0.5*b+rec(1)+(j-1)*bp/ndt+0.5*(bp/ndt);
        eMecProfile(i,2)=x;
        
        y=-0.5*h+rec(2)+0.5*t2;
        eMecProfile(i,3)=y;
    elseif i>ndt*2 && i<=ndt*3
        j=i-2*ndt; %counting the number of dt along a face
        eMecProfile(i,1)=dA3;
        x=-0.5*b+rec(1)+0.5*t3;
        eMecProfile(i,2)=x;
        
        y=0.5*h-rec(2)-(j-1)*hp/ndt-0.5*(hp/ndt);
        eMecProfile(i,3)=y;
    elseif i>ndt*3 && i<=ndt*4
        j=i-3*ndt;
        eMecProfile(i,1)=dA4;
        x=0.5*b-rec(1)-0.5*t4;
        eMecProfile(i,2)=x;
        
        y=0.5*h-rec(2)-(j-1)*hp/ndt-0.5*(hp/ndt);
        eMecProfile(i,3)=y;
    end
end

nv=4*ndt;
sumaM=0;
sumaF=0;
for i=1:nv
    eMecProfile(i,4)=0.5*h-eMecProfile(i,3); %momentum distance
    eMecProfile(i,5)=0.003/c*(eMecProfile(i,4)-c);
    
    if (eMecProfile(i,5)<-0.0021)
        eMecProfile(i,5)=-0.0021;
    elseif(eMecProfile(i,5)>0.0021)
        eMecProfile(i,5)=0.0021;
    end
    eMecProfile(i,6)=eMecProfile(i,5)*E;
    eMecProfile(i,7)=eMecProfile(i,6)*eMecProfile(i,1);
    eMecProfile(i,8)=eMecProfile(i,7)*(eMecProfile(i,4)-0.5*h);
    
    sumaF=sumaF+eMecProfile(i,7);
    sumaM=sumaM+eMecProfile(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcreto(a,fdpc,b,h);

eleMec=[elemAc;
        elemConc];
      