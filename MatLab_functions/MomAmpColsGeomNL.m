function [Delta,Mamp]=MomAmpColsGeomNL(fc,k,I,L,V,P,M,b,h,plotdef)

%------------------------------------------------------------------------
% Syntax:
% [Delta,Pcr,Mamp]=MomAmpColsGeomNL(fc,k,I,L,V,P,Md,b,h,plotdef)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To compute the amplified moments for a column considering the
%          P-Delta effects by geometric Non-Linearity
% 
% OUTPUT: Delta:    Lateral displacement or additional load eccentricity
%                   that causes the second-order moments
%
%         Mamp:     amplified moment in the direction of question
%
%         Pcr:      is the critical axial load of Euler for instability
%
% INPUT:  P:        Is the axial load over the column's cross-section
%
%         V:        Is the shear force at the top of the column
%
%         fc:       is the f'c used
%
%         b,h:      are the column cross-section dimensions
%
%         M:        is the acting bending moment:
%           
%         k:        is the slenderness factor, according to the boundary
%                   conditions
%
%         I:        Momentum of inertia of the cross-section in the current
%                   axis of reference
%
%         plotdef:  is the parameter that indicates if the plot of the 
%                   deformed column is required
%
%         L:        total length of the column element (without 
%                   considering the slenderness factor)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-07-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% Topology -----------------------------------------------------------
nnodes=2;
Edof=[1  1  2  3  4  5  6]; 

%% Element properties and global coordinates -------------------------- 
if fc<2000 % units: (kg,cm)
    E=14000*sqrt(fc);  
else       % units: (lb,in)
    E=57000*sqrt(fc);  
end

A1=b*h;       
I1=I;   
ep1=[E A1 I1];	 
eq1=[0];       
ex1=[0 0];         
ey1=[0 L];     

%% Initial axial forces ----------------------------------------------

QX1=0.0001; 
QX01=1;

%% Boundary conditions -----------------------------------------------

% Restricted DODF in function of the slenderness factor
if k==2 
    bc=[1 0;2 0;3 0];	
elseif k==0.5
    bc=[1 0;2 0;3 0;6 0];
elseif k==1
    bc=[1 0;2 0;4 0];
end

%% Iteration for convergence -----------------------------------------
ecc=abs(M/P);
eps=1e-6;		
n=0;			
while(abs((QX1-QX01)/QX01)>eps)
    n=n+1;

    K=zeros(3*nnodes,3*nnodes);
    f=zeros(3*nnodes,1);	
    f(4)=V; f(5)=P; f(6)=abs(P*ecc);

    [Ke1,fe1]=beam2ge(ex1,ey1,ep1,QX1,eq1); % element stiffness matrix

    [K,f]=assem(Edof(1,:),K,Ke1,f,fe1); % assembles stiffness matrix

    [displacements_Nonlinear,r]=solveq(K,f,bc);
    EdNonlin=extract_ed(Edof,displacements_Nonlinear);
    QX01=QX1; 
    [esNonlinear,QX1,ediNonLin]=beam2gs(ex1,ey1,ep1,EdNonlin(1,:),QX1,eq1,7);
    if(n==1)
        edilin=ediNonLin;
        eslinear=esNonlinear;
        Edlinear=EdNonlin;
        displacements_linear=displacements_Nonlinear;
    end  

    if(n>20)
        disp('The solution does not converge')
        break
    end
end
Mamp=max(abs(esNonlinear(:,3)));
Delta=Mamp/abs(P)-ecc;

%% Draw deformed frame ---------------------------------------------
if plotdef==1
    figure(1)
    plotpar=[3 1 0];
    eldraw2(ex1,ey1,plotpar);
    sfac=scalfact2(ex1,ey1,ediNonLin,0.1);
    plotpar=[1 2 0];
    dispbeam2(ex1,ey1,ediNonLin,plotpar,sfac);
    plotpar=[2 4 0];
    dispbeam2(ex1,ey1,edilin,plotpar,sfac);
    title('Displacements')
end

%------------------------ end --------------------------------------
