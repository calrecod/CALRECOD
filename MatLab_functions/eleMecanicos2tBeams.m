function eleMec=eleMecanicos2tBeams(c,a,fdpc,h,b,b_rec,h_rec,E,t1,t2)

%------------------------------------------------------------------------
% Syntax:
% eleMec=eleMecanicos2tBeams(c,a,fdpc,h,b,b_rec,h_rec,E,t1,t2)
%
%------------------------------------------------------------------------
% PURPOSE: To compute the sum of resistant forces of a beam cross-section,
% considering the contribution of steel in tension, steel in compression 
% and concrete in compression.
% 
% OUTPUT: eleMec:       vector that contains the output [Fs, Ms;
%                                                        Fc, Mc]
%
% INPUT:  t1,t2         are the given width of ISR in compression and 
%                       tension, respectively
%
%         b_rec,h_rec:  are the concrete cover parameters horizontally and
%                       vertically, respectively
%
%         fdpc:         Is the reduced f'c as 0.85f'c according to ACI 318
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bp=b-2*b_rec;

eMecProfile=zeros(2,8);

% eMecProfile=[dA, xposition, yposition, d, eps, eE, Fr, Mr]

%calculate dA for each t area
A1=(b-2*b_rec)*t1;
A2=(b-2*b_rec)*t2;

% localtion of each discrete dt over the section
for i=1:2
    if i==1
        j=i; %counting the number of dt along a face
        eMecProfile(i,1)=A1;
        x=-0.5*b+b_rec+(j-1)*bp+0.5*(bp);
        eMecProfile(i,2)=x;
        
        y=0.5*h-h_rec-0.5*t1;
        eMecProfile(i,3)=y;
        
    elseif i==2
        j=i-1; %counting the number of dt along a face
        eMecProfile(i,1)=A2;
        x=-0.5*b+b_rec+(j-1)*bp+0.5*(bp);
        eMecProfile(i,2)=x;
        
        y=-0.5*h+h_rec+0.5*t2;
        eMecProfile(i,3)=y;
    end
end

nv=2;
sumaM=0;
sumaF=0;
for i=1:nv
    eMecProfile(i,4)=0.5*h-eMecProfile(i,3); %d1
    eMecProfile(i,5)=0.003/c*(eMecProfile(i,4)-c); %epsilum
    
    if(eMecProfile(i,5)>0.0021)
        eMecProfile(i,5)=0.0021; % Only the steel in compression
                                 % is limited by its yielding strain
    end
    eMecProfile(i,6)=eMecProfile(i,5)*E;
    eMecProfile(i,7)=eMecProfile(i,6)*eMecProfile(i,1)*0.001;
    eMecProfile(i,8)=eMecProfile(i,7)*(eMecProfile(i,4)-0.5*h)*0.01;
    
    sumaF=sumaF+eMecProfile(i,7);
    sumaM=sumaM+eMecProfile(i,8);
end

elemAc=[sumaF sumaM];
elemConc=casoConcreto(a,fdpc,b,h);

eleMec=[elemAc;
        elemConc];
end
      