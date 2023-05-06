function [bestPerformance,bestPosition,besttippingFS,bestslideFS,...
    bestLCap_FS,bestsepheel,bestefHeel,bestsepfoot,besteffoot,bestseptrunk,...
    besteftrunk]=DesignRetainingWallPSO(minDim,MaxDim,H,D,m1,m2,FiFill,...
    wvFill,beta,FiBackFill,alfa,FiFound,fc,fy,wvc,ductility,qadm,minFSqadm,...
    SlideSF,TippingSF,typeRebar,sepMinRebars,maxEf,qaf,qab,qs,LF_DL)
        
%-------------------------------------------------------------------------
% Syntax:
% [bestPerformance,bestPosition,besttippingFS,bestslideFS,...
%  bestLCap_FS,bestsepheel,bestefHeel,bestsepfoot,besteffoot,bestseptrunk,...
%  besteftrunk]=DesignRetainingWallPSO(minDim,MaxDim,H,D,m1,m2,FiFill,...
%  wvFill,beta,FiBackFill,alfa,FiFound,fc,fy,wvc,ductility,qadm,minFSqadm,...
%  SlideSF,TippingSF,typeRebar,sepMinRebars,maxEf,qaf,qab,qs,LF_DL)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%
%-------------------------------------------------------------------------
% PURPOSE: To optimally design a reinforced concrete retaining wall with
% the PSO algorithm. The design variables "foot", "heel", "hf" and "b" are
% taken as the optimization variables. The variables m1 and m2 
% (corresponding to the front and back wall's trunk slopes) remain
% constant thorughout the process.
% 
% OUTPUT: bestPerformance:      is the wall's final cross-sectional area of
%                               the optimal design
%
%         besttippingFS:        is the final tipping Safety Factor for the
%                               optimally designed wall
%
%         bestslideFS:          is the final slide Safety Factor for the
%                               optimally designed wall
%
%         bestLCap_FS:          is the final Safety Factor against the
%                               soil's bearing load capacity for the
%                               optimally designed wall
%
%         bestsepheel:          is the final rebar separation for the
%                               designed reinforcement in the optimally 
%                               designed concrete wall's heel
%
%         bestefHeel:           is the final structural efficiency for the
%                               optimally designed wall's heel
%
%         bestsepfoot:          is the final rebar separation for the
%                               designed reinforcement in the optimally 
%                               designed concrete wall's foot
%
%         besteffoot:           is the final structural efficiency for the
%                               optimally designed wall's foot
%
%         bestseptrunk:         is the final rebar separation for the
%                               designed reinforcement in the optimally 
%                               designed concrete wall's trunk
%
%         besteftrunk:          is the final structural efficiency for the
%                               optimally designed wall's trunk
%
% INPUT:  minDim:               is the vector containing the min dimension
%                               values that the optimization variables can 
%                               have during the optimization process, in 
%                               format: [foot,heel,hf,b]
%
%         MaxDim:               is the vector containing the max dimension
%                               values that the optimization variables can 
%                               have during the optimizatoin process, in 
%                               format: [foot,heel,hf,b]
%
%         FiFill:               is the fill soil's friction angle
%
%         H:                    is the wall's stem height dimension
%
%         D:                    is the soil's bakc fill depth
%
%         m1, m2:               are the wall's dowels' from and back slopes
%
%         wvFill:               is the soil fill's volumetric weight
%
%         beta:                 is the front fill soil's upper-grade angle
%
%         FiBackFill:           is the back fill soil's friction angle
%
%         alfa:                 is the back fill soil's upper-grade angle
%
%         FiFound:              is the foundation soil's friction angle
%
%         fc:                   is the concrete compressive strength f'c
%                               (in units Kg/cm^2)
%
%         wvc:                  is the concrete volumetric weight (Kg/cm^3)
%
%         ductility:            is the level of ductility demand for the
%                               design of the reinforcement
%
%         qadm:                 is the soil's bearing load capacity, in 
%                               units (Kg/cm^2)
%
%         minFSqadm:            is the design Safety Factor against the
%                               soil's bearing load capacity
%
%         SlideSF:              is the design Safety Factor against the
%                               slide forces over the wall
%
%         TippingSF:            is the design Safety Factor against the
%                               tipping forces over the wall
% 
%         typeRebar:            is the rebar eigth-of-an-inch to be used
%                               for the wall's reinforcement
%
%         sepMinRebars:         is the minimum rebar separation restriction
%                               for the wall's reinforcement
%
%         maxEf:                is the critical structural efficiency for 
%                               the design of the wall's heel, foot and 
%                               trunk
%
%         qaf, qab:             is the front surcharge and the back 
%                               surcharge, respectively (in units Kg/cm^2) 
%
%         qs:                   is the linear load that the wall may
%                               support over the top along its length
%
%         LF_DL:                is the Dead Load Design Factor
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-04-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% PSO algorithmic parameters 
alphaPSO=1;
c1=2; % cognitive component
c2=2; % social component
dt=0.5;
inertiaWeight=1.3;
betaPSO=0.99;

numberOfParticles=50;
numberOfDimensionSpace=4;
niter=40;

%% Optimization variables 

maxVelocity=(MaxDim-minDim)./dt;

%% Generate position and velocity vector of each particle
    
for i=1:numberOfParticles
    for j=1:numberOfDimensionSpace
        r=rand;
        PositionMatrix(i,j)=minDim(j)+r*(MaxDim(j)-minDim(j));
        velocityMatrix(i,j)=alphaPSO/dt*(-(MaxDim(j)-minDim(j))*0.5+r*(MaxDim(j)-minDim(j)));
    end
end


%% Optimization process
bestPerformance=inf;
bestPositionSwarmMatrix=PositionMatrix;
bestPerformanceSwarm=inf(numberOfParticles,1);
for iter=1:niter
    
    % Determine the best position and best performance
    for i=1:numberOfParticles
        position=PositionMatrix(i,:);
        foot=position(1);
        heel=position(2);
        hf=position(3);
        b=position(4);
        
        %% Objective function
        [compliedRestric,performance(i),linearWeigthWall,tippingFS,slideFS,...
        LCap_FS,sepheel,efHeel,sepfoot,effoot,septrunk,eftrunk]=RetainingRCWall...
        (foot,heel,hf,b,FiFill,H,D,m1,m2,wvFill,beta,FiBackFill,alfa,FiFound,...
        fc,fy,wvc,ductility,qadm,minFSqadm,SlideSF,TippingSF,typeRebar,...
        sepMinRebars,maxEf,qaf,qab,qs,LF_DL);
    
        if (performance(i)<bestPerformance && compliedRestric==1)
            bestPerformance=performance(i);
            bestParticlePerformanceIndex=i; %to store later the globalbest
            bestPosition=PositionMatrix(i,:);
            
            besttippingFS=tippingFS;
            bestslideFS=slideFS;
            bestLCap_FS=LCap_FS;
            bestsepheel=sepheel;
            bestefHeel=efHeel;
            bestsepfoot=sepfoot;
            besteffoot=effoot;
            bestseptrunk=septrunk;
            besteftrunk=eftrunk;
        end
        
        if (performance(i)<bestPerformanceSwarm(i) && compliedRestric==1)
            bestPerformanceSwarm(i)=performance(i);
            bestPositionSwarmMatrix(i,:)=position;
        end
    end
    
    %% Global best position 
    globalBestPosition=bestPosition;
    globalBestPerformance=bestPerformance;
    
    %% Velocity and position uptdating 
    
    for i=1:numberOfParticles
        q=rand;
        r=rand;
        
        % Velocity updating
        for j=1:numberOfDimensionSpace
            velocityMatrix(i,j)=inertiaWeight*velocityMatrix(i,j)+...
                c1*q*((bestPositionSwarmMatrix(i,j)-PositionMatrix(i,j))/dt)+...
                c2*r*((bestPosition(j)-PositionMatrix(i,j)));
        end
        absolouteValueVelocity=abs(velocityMatrix(i,:));
        valueVelocity=velocityMatrix(i,:);
        
        for j=1:numberOfDimensionSpace
            if (absolouteValueVelocity(j)>maxVelocity(j))
                    velocityMatrix(i,j)=maxVelocity(j);

            end
            if (valueVelocity(j)<-maxVelocity(j))
                    velocityMatrix(i,j)=-maxVelocity(j);

            end
            
        end

        % position uptdating
        for j=1:numberOfDimensionSpace
            newPosition=abs(PositionMatrix(i,j)+...
                                velocityMatrix(i,j)*dt);
            if newPosition<minDim(j)
                newPosition=minDim(j);
            elseif newPosition>MaxDim(j)
                newPosition=MaxDim(j);
            end
            newPosition=newPosition-mod(newPosition,5)+5;
            PositionMatrix(i,j)=newPosition;
        end
    end
    if(inertiaWeight*betaPSO<0.3)
        inertiaWeight=inertiaWeight;
    else
        inertiaWeight=inertiaWeight*betaPSO;
    end

end
