function [da4,coordISRdisc]=discreteISRColRec(t,b,h,rec,ndt)
%------------------------------------------------------------------------
% Syntax:
% [da4,coordISRdisc]=discreteISRColRec(t,b,h,rec,ndt)
%------------------------------------------------------------------------
% SYSTEM OF UNITS: Any
% 
%------------------------------------------------------------------------
% PURPOSE: To discretize the ISR of a rectangular column cross-section
% according to the number of discrete elements ndt.
%
% Note: each of the four ISR cross-section's boundaries will be discretized
% in ndt elements.
% 
% OUTPUT: coordISRdisc:        Are the local coordinates of the discrete
%                              ISR elements over the rectangular cross-
%                              section
%
%         da4:                 is the cross-section area of each discrete
%                              ISR's element
% 
% INPUT:  b,h:                cross-section dimensions of column (width 
%                             and height)
%
%         t:                  is the ISR width of a 1t-ISR for column 
%                             cross-sections
%
%         rec:                is a vector containing the concrete cover for
%                             both cross-section axis as [coverX,coverY]
%
%         ndt:                is the number of elements in which the ISR is
%                             to be divided for each of its four boundaries
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-22
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

bp=b-2*rec(1);
hp=h-2*rec(2);

% calculate dA for each t area
dA1=(bp)*t/ndt;
dA2=(bp)*t/ndt;
dA3=(hp)*t/ndt;
dA4=(hp)*t/ndt;

da4=[dA1,dA2,dA3,dA4];

% localtion of each discrete dt over the section
for i=1:ndt*4
    if i<=ndt*1
        j=i; % counting the number of dt along the upper cross-section 
             % boundary
        x=-0.5*b+rec(1)+(j-1)*bp/ndt+0.5*(bp/ndt);
        coordISRdisc(i,1)=x;
        
        y=0.5*h-rec(2)-0.5*t;
        coordISRdisc(i,2)=y;
        
    elseif i>ndt*1 && i<=ndt*2
        j=i-ndt; % counting the number of dt along the upper cross-section 
                 % boundary
        x=-0.5*b+rec(1)+(j-1)*bp/ndt+0.5*(bp/ndt);
        coordISRdisc(i,1)=x;
        
        y=-0.5*h+rec(2)+0.5*t;
        coordISRdisc(i,2)=y;
    elseif i>ndt*2 && i<=ndt*3
        j=i-2*ndt; % counting the number of dt along the upper 
                   % cross-section boundary
        x=-0.5*b+rec(1)+0.5*t;
        coordISRdisc(i,1)=x;
        
        y=0.5*h-rec(2)-(j-1)*hp/ndt-0.5*(hp/ndt);
        coordISRdisc(i,2)=y;
    elseif i>ndt*3 && i<=ndt*4
        j=i-3*ndt; % counting the number of dt along the upper 
                   % cross-section boundary
        x=0.5*b-rec(1)-0.5*t;
        coordISRdisc(i,1)=x;
        
        y=0.5*h-rec(2)-(j-1)*hp/ndt-0.5*(hp/ndt);
        coordISRdisc(i,2)=y;
    end
end