function [Mr_col,h,Inertia_xy_modif,bestPerformance,best_cost,bestnv,...
    bestPerformanceEf,bestArrangement,best_disposicion]=PSOAsymmetricRebar...
    (t4,b,h,rec,Ac_t_elem,npuntos,fdpc,beta1,pu_asym_cols,load_conditions,...
    condition_cracking,plotRebarDesign,rebarOptimConv)

%------------------------------------------------------------------------
% Syntax:
% [Mr_col,h,Inertia_xy_modif,bestPerformance,best_cost,bestnv,...
%   bestPerformanceEf,bestArrangement,best_disposicion]=PSOAsymmetricRebar...
%   (t4,b,h,rec,Ac_t_elem,npuntos,fdpc,beta1,pu_asym_cols,load_conditions,...
%   condition_cracking,plotRebarDesign,rebarOptimConv)
%
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal arrangement of rebars asymmetrically 
% over a column cross-section using the PSO algorithm based on a given 
% 1t-ISR.
% 
% OUTPUT: Mr_col:                   are the resisting bending moments for 
%                                   both axis directions for the optimally 
%                                   desined cross-section: [M{Rx},M{Ry}]
%
%         h:                        is the increased height dimension of 
%                                   the cross-section; in case it suffers
%                                   modification through the optimization 
%                                   design process to comply the 
%                                   corresponding restrictions
%
%         Inertia_xy_modif:         are the modified inertia momentums for
%                                   both axis directions of the optimally
%                                   reinforced cross-section: under cracking
%                                   criteria (non-cracked or cracked)
%
%         bestPerformance:          is the optimal rebar area for the 
%                                   cross-section
%
%         best_cost:                is the final construction cost for the
%                                   optimally designed cross-section 
%                                   (considering only rebar volumes and 
%                                   assembly)
%
%         bestnv:                   is the total number of rebars to be 
%                                   placed over the optimally designed 
%                                   cross-section
%
%         bestPerformanceEf:        is the critical structural efficiency
%                                   for the optimally designed cross-section
%
%         bestArrangement:          is the list of rebar types for each of
%                                   the rebars to be placed over the 
%                                   optimally designed cross-section
%
%         best_disposicion:         are the local coordinates of the optimal
%                                   rebar option
%
% INPUT:  fdpc:                 is the f'c reduced with the factor 0.85 
%                               according to code
%
%         beta1:                is determined as stablished by code (see
%                               Documentation)
%
%         t4:                   is the optimal ISR from which the optimal 
%                               rebar option is determined as: [t1,t2,t3,t4]
%                               (see Dcoumentation)
%
%         b,h:                  cross-section initial dimensions
%
%         npuntos:              number of points to be computed for the
%                               definition of the interaction diagram
%
%         sepMin:               is the minimum rebar separation constriction 
%                               (see Documentation)
%
%         rec:                  are the concrete cover values for each axis
%                               direction of the cross-section
%
%         pu_asym_cols:         are the unit cost data for the reinforcing 
%                               bar volumes and assembly as an array of size:
%                               [2,7] by default, for which the first row
%                               corresponds to unit-cost values for each 
%                               rebar type in case only one type of rebar
%                               results as an optimal design, and the second 
%                               row consisting of only one unit cost in case 
%                               more than one different type of rebar results
%                               as an optimal design (assuming an average
%                               unit-cost of all types of rebars)
%
%         load_conditions:      is the load condition array: size=[nloads,4]
%                               in format [nload,Pu,Mux,Muy]
%
%         condition_cracking:   is the parameter to indicate the cracking
%                               mechanism to be considered; options are 
%                               ''Cracked'' or ''Non-cracked''
%
%         plotRebarDesign:      is the parameters that indicates if the 
%                               rebar design results are required or not. 
%                               Options are: (1) they are required, (2) 
%                               they are not required
%
%         rebarOptimConv:       is the parameters that indicates if the 
%                               optima rebar convergence is required or not.
%                               Options are: (1) they are required, (2) 
%                               they are not required
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

Mr_col=zeros(1,2);
fy=4200;
E=2e6;
nopciones=0;
while nopciones==0
    
    bestInteracDiagram=zeros(npuntos+1,9);
    bp=b-2*rec(1);
    hp=h-2*rec(2);

    maxEfficency=1.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Particle Swarm Optimization Algorithm parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha=1;
    c1=2; % cognitive component
    c2=2; % social component
    dt=1.0;
    inertiaWeight=1.3;
    beta=0.99;

    %%%% Posible reinforcement combos _________________________________

    %__ Data base rebar ...............................................
    varDisponibles=[4 4/8*2.54 0.994;
                    5 5/8*2.54 1.552;
                    6 6/8*2.54 2.235;
                    8 8/8*2.54 3.973;
                    9 9/8*2.54 5.033;
                    10 10/8*2.54 6.225;
                    12 12/8*2.54 8.938];

    a1=t4(1)*(b-2*rec(1));
    a2=a1;
    a3=t4(3)*(h-2*rec(2));
    a4=a3;

    number_rebars_sup=zeros(1,7);
    number_rebars_inf=zeros(1,7);
    number_rebars_izq=zeros(1,7);
    number_rebars_der=zeros(1,7);

    for i=1:7
        if fix(a1/(pi/4*(varDisponibles(i,2))^2))+1<2
            number_rebars_sup(i)=2;
        else
            number_rebars_sup(i)=fix(a1/(pi/4*(varDisponibles(i,2))^2))+1;
        end

        if fix(a2/(pi/4*(varDisponibles(i,2))^2))+1<2
            number_rebars_inf(i)=2;
        else
            number_rebars_inf(i)=fix(a2/(pi/4*(varDisponibles(i,2))^2))+1;
        end

        number_rebars_der(i)=fix(a4/(pi/4*(varDisponibles(i,2))^2))+1;
        number_rebars_izq(i)=fix(a3/(pi/4*(varDisponibles(i,2))^2))+1;

    end
    
    
    %_____________________________________________________________________
    numberOfParticles=25;
    numberOfDimensionSpace=4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%Generate position and velocity vector of each particle%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ats=zeros(numberOfParticles,1);
    globalBestAreaMatrix=zeros(numberOfParticles,1);

    PositionMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
    velocityMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
    for i=1:numberOfParticles

        for j=1:numberOfDimensionSpace

            r=rand;
            xmax=7;
            xmin=1;
            PositionMatrix(i,j)=fix(xmin+r*(xmax-xmin))+1;

            velocityMatrix(i,j)=alpha/dt*(-(xmax-xmin)*0.5+...
                r*(xmax-xmin));
        end
    end

    %_______________________________________________________________
    %Main Loop
    %_______________________________________________________________

    iterationsNumber=35;

    % Se busca maximizar la eficiencia (máx 0.95) o 95 por ciento
    bestPerformance=0;
    bestPerformanceEf=0;
    minAts=inf;

    bestPosition=zeros(1,numberOfDimensionSpace);
    position=zeros(1,numberOfDimensionSpace);
    performance=zeros(numberOfParticles,1);
    bestPositionSwarmMatrix=PositionMatrix;
    bestPerformanceSwarm=zeros(numberOfParticles,1);
    GlobalBestEfPerformanceMatrix=inf(numberOfParticles,1);
    globalBestPerformance=0;
    globalBestAreaMatrix=zeros(numberOfParticles,1);
    for iteration=1:iterationsNumber

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the best position and best performance %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for i=1:numberOfParticles
            position=PositionMatrix(i,:);

            op1=position(1);
            op2=position(2);
            op3=position(3);
            op4=position(4);
            dbmin=max([varDisponibles(op1,2),varDisponibles(op2,2),...
                      varDisponibles(op3,2),varDisponibles(op4,2)]);

            sepMin=max([1.5*2.54 3/2*dbmin]);

            Ats(i)=number_rebars_sup(op1)*(varDisponibles(op1,2)^2*pi/4)+...
                    number_rebars_inf(op2)*(varDisponibles(op2,2)^2*pi/4)+...
                    number_rebars_izq(op3)*(varDisponibles(op3,2)^2*pi/4)+...
                    number_rebars_der(op4)*(varDisponibles(op4,2)^2*pi/4);

            At=Ats(i);
            nv=number_rebars_sup(op1)+number_rebars_inf(op2)+...
                number_rebars_izq(op3)+number_rebars_der(op4);
            
            if op1==op2 && op2==op3 && op3==op4
                cost_Rebar_Asymm=At*0.0001*7800*(pu_asym_cols(1,op1));
            else
                cost_Rebar_Asymm=At*0.0001*7800*pu_asym_cols(2,1);
            end
            
            [disposition_rebar,separacion_hor1,separacion_hor2,...
            separacion_ver1,separacion_ver2]=dispositionRebarAsymmetric(b,...
            h,sepMin,rec,nv,number_rebars_sup,number_rebars_inf,...
            number_rebars_izq,number_rebars_der,varDisponibles,op1,op2,op3,op4);
            
            [performance(i),diagramaInteraccion,eficiencia_table,cp,cxy]=...
                EvaluateAsymmetric(load_conditions,npuntos,position,b,h,...
                fy,fdpc,beta1,E,number_rebars_sup,number_rebars_inf,...
                number_rebars_izq,number_rebars_der,varDisponibles,disposition_rebar);

            if (Ats(i)<=minAts && performance(i)<maxEfficency && ...
                    separacion_hor1>=sepMin && separacion_hor2>=sepMin && ...
                    separacion_ver1>=sepMin && separacion_ver2>=sepMin)
                
                bestSepRebars=[separacion_hor1,separacion_hor2,...
                                separacion_ver1,separacion_ver2];
                best_cost=cost_Rebar_Asymm;
                bestMrx=eficiencia_table(1,5);
                bestMry=eficiencia_table(1,7);
                minAts=Ats(i);
                best_cp=cp;
                best_cxy=cxy;
                bestPerformance=Ats(i);
                bestInteracDiagram=diagramaInteraccion;
                bestPerformanceEf=performance(i);
                bestParticlePerformanceIndex=i; %to store later the globalbest
                bestPosition=position;
                bestnv=nv;
                best_disposicion=disposition_rebar;
            end   

            % global best ------------
            % here past performance is compared to deterine
            % a global best entity
            if (Ats(i)<=globalBestAreaMatrix(i) && ...
                    performance(i)<maxEfficency && ...
                    separacion_hor1>=sepMin && separacion_hor2>=sepMin && ...
                    separacion_ver1>=sepMin && separacion_ver2>=sepMin)

                bestPerformanceSwarm(i)=Ats(i);
                bestPositionSwarmMatrix(i,:)=position;
                globalBestAreaMatrix(i)=Ats(i);

                GlobalBestEfPerformanceMatrix(i)=performance(i);
            end
           globalBestPosition=bestPosition;
           globalBestPerformance=bestPerformance;

        end
        
        % PLOTTING EVOLUTION OF OPTIMIZATION CONVERGENCE________________
        if rebarOptimConv==1
            figure(10)
            plot(iteration,...
                bestPerformance,'b o','linewidth',0.1,'MarkerFaceColor','blue' )
            hold on
            xlabel('Iteration')
            ylabel('Area (cm2)')
            title('Global Optimal Reinforcing Bar Area')
            hold on
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % velocity and position uptdating %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for i=1:numberOfParticles
            q=rand;
            r=rand;
            % velocity updating
            for j=1:numberOfDimensionSpace
                velocityMatrix(i,j)=inertiaWeight*velocityMatrix(i,j)+...
                    c1*q*((bestPositionSwarmMatrix(i,j)-PositionMatrix(i,j))/dt)+...
                    c2*r*((bestPosition(j)-PositionMatrix(i,j)));

                xmax=7;
                xmin=1;
                maxVelocity=(xmax-xmin)/dt;
                absolouteValueVelocityX=abs(velocityMatrix(i,j));
                valueVelocityX=velocityMatrix(i,j);

                if (absolouteValueVelocityX>maxVelocity)
                    velocityMatrix(i,j)=maxVelocity;
                end
                if (valueVelocityX<-maxVelocity)
                    velocityMatrix(i,j)=-maxVelocity;
                end
            end


            %%%%%%%%%%%%%%%%%%%%
            % position updating% 
            % t width updating %
            %%%%%%%%%%%%%%%%%%%%
            Ats(i)=0;
            for j=1:numberOfDimensionSpace
                PositionMatrix(i,j)=fix(PositionMatrix(i,j)+velocityMatrix(i,j)*dt);
                if PositionMatrix(i,j)<=0
                    PositionMatrix(i,j)=1;
                elseif PositionMatrix(i,j)>7
                    PositionMatrix(i,j)=7;
                end
            end
        end

        if(inertiaWeight*beta<0.3)
            inertiaWeight=inertiaWeight;
        else
            inertiaWeight=inertiaWeight*beta;
        end
    end

    
    if sum(bestPosition)==0
        nopciones=0;
        fprintf('\nDimensiones no aptas para columna con refuerzo\n');
        fprintf('asimétrico; Se aumentará el peralte 5cm\n');
        h=h+5;
    else 
        nopciones=1;
    end
   
end

bestop1=bestPosition(1);
bestop2=bestPosition(2);
bestop3=bestPosition(3);
bestop4=bestPosition(4);

bestnv1=number_rebars_sup(bestop1);
bestnv2=number_rebars_inf(bestop2);

bestnv3=number_rebars_izq(bestop3);
bestnv4=number_rebars_der(bestop4);

%fprintf('Best option Rebar combo \n %.6f\t%.6f\t%.6f\t%.6f\n',bestPosition);
%fprintf('\nBest area=%.6f\n',bestPerformance);
%fprintf('Best efficiency=%.2f\n',bestPerformanceEf);
%fprintf('Best cost=%.2f\n',best_cost);
%fprintf('Best CP=%.2f\n',best_cp);

Mr_col(1)=bestMrx;
Mr_col(2)=bestMry;

pu=load_conditions(1,2);
mux=load_conditions(1,3);
muy=load_conditions(1,4);

excentricity_x=abs(mux/pu)*100;
excentricity_y=abs(muy/pu)*100;
eccentricity_xy=[excentricity_x,excentricity_y];

t1_var=bestnv1*(varDisponibles(bestop1,2)^2*pi/4)/(b-2*rec(1)); 
t2_var=bestnv2*(varDisponibles(bestop2,2)^2*pi/4)/(b-2*rec(1));
t3_var=bestnv3*(varDisponibles(bestop3,2)^2*pi/4)/(h-2*rec(2));
t4_var=bestnv4*(varDisponibles(bestop4,2)^2*pi/4)/(h-2*rec(2));

%%% Modified inertia momentum for the reinforced cross-seciton............
[Inertia_xy_modif,Atransf_xy]=CrackingColumnsAsym(h,b,fdpc,rec,eccentricity_xy,...
                            t1_var,t2_var,t3_var,t4_var,pu,best_cxy,condition_cracking,best_cp);
              
bestArrangement=zeros(bestnv,1);
bestArrangement(1:bestnv1)=bestop1;
bestArrangement((bestnv1+1):(bestnv1+bestnv2))=bestop2;
bestArrangement((bestnv1+bestnv2+1):(bestnv1+bestnv2+bestnv3))=bestop3;
bestArrangement((bestnv1+bestnv2+bestnv3+1):(bestnv1+bestnv2+bestnv3+bestnv4))=bestop4;

if plotRebarDesign==1
    diagramsFinalRebarCols(load_conditions,bestInteracDiagram,best_disposicion,...
                        h,b,bestArrangement);
end

