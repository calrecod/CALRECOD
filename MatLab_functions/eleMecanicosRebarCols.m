function [eMecVar]=eleMecanicosRebarCols(disposition_rebar,nv,ov,av,b,h,...
                                      c,fdpc,E,beta1)

%------------------------------------------------------------------------
% Syntax:
% [eMecVar]=eleMecanicosRebarCols(disposition_rebar,nv,ov,av,b,h,...
%                                     c,fdpc,E,beta1)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: To compute the sum of resistant forces of a reinforced column 
% cross-section considering the distribution of rebars over the
% cross-section and concrete zone in compression.
% 
% OUTPUT: eMecVar:              vector that contains the output [Fs,Ms;
%                                                                Fc,Mc]
%
% INPUT:  c:                    is the neutral axis depth
%
%         E:                    Elasticity modulus of steel
%
%         b,h:                  cross-section dimensions
%
%         fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta1:                is determined as stablished by code (see
%                               Documentation)
%
%         nv:                   is the number of rebars to be placed over
%                               the cross-section
%
%         ov,av:                are the type of rebar in eighth of inches 
%                               (ov/8 in) and the cross-section area of each
%                               rebar
%
%         rebar_disposition:    are the local coordinates of rebars over 
%                               the cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

eMecVarillado=zeros(nv,9);
eMecVarillado(:,1)=ov;
eMecVarillado(:,2)=disposition_rebar(:,1);
eMecVarillado(:,3)=disposition_rebar(:,2);

sumaM=0;
sumaF=0;
        
for i=1:nv
    eMecVarillado(i,4)=0.5*h-eMecVarillado(i,3);
    eMecVarillado(i,5)=av;
    eMecVarillado(i,6)=0.003/c*(eMecVarillado(i,4)-c);
    if (eMecVarillado(i,6)<-0.0021)
        eMecVarillado(i,6)=-0.0021;
    elseif(eMecVarillado(i,6)>0.0021)
        eMecVarillado(i,6)=0.0021;
    end
    eMecVarillado(i,7)=eMecVarillado(i,6)*E;
    eMecVarillado(i,8)=eMecVarillado(i,7)*av;
    eMecVarillado(i,9)=eMecVarillado(i,8)*(eMecVarillado(i,4)-0.5*h);
    
    sumaF=sumaF+eMecVarillado(i,8);
    sumaM=sumaM+eMecVarillado(i,9);
end

a=beta1*c;
elemConc=casoConcreto(a,fdpc,b,h);
eMecVar=[sumaF sumaM;
    elemConc(1) elemConc(2)];

end