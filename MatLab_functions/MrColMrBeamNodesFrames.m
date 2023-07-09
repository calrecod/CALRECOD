function [EffMcb]=MrColMrBeamNodesFrames(Mp,elemcols,elembeams,ni,nf,nnodes)
%-------------------------------------------------------------------------
% SYNTAX : 
% [EffMcb]=MrColMrBeamNodesFrames(Mp,elemcols,elembeams,ni,nf,nnodes)
%-------------------------------------------------------------------------
%    PURPOSE
%     To revise the criteria "Strong Column - Weak Beam" of each node of 
%     a 2D/3D frame by summing up the intersecting beams' resistant moment
%     and the intersecting columns' resistant moment.
% 
%    OUTPUT:  EffMcb:                   is the beam-column node bending
%                                       efficiency: SumMrCol / SumMrBeam
%
%    INPUT:   ni,nf:                    vector containing the initial and
%                                       final nodes of all elements.
%                                       Size: [n-elems, 1]
% 
%             elembeams:                vector containing the list of 
%                                       elements classified as beams
% 
%             elemcols:                 vector containing the list of 
%                                       elements classified as columns
%
%             nnodes:                   number of nodes of the 2D Frame
%
%-------------------------------------------------------------------------
% LAST MODIFIED: L.Verduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%-------------------------------------------------------------------------
beams=length(elembeams);
cols=length(elemcols);

EffMcb=zeros(nnodes,1);
MbNode=zeros(nnodes,1);
McNode=zeros(nnodes,1);
for node=1:nnodes
    % To detect which beam element intersect with the node (if any)
    sumMrBeams=0;
    for j=1:beams
        
        elb=elembeams(j);
        if ni(elb)==node
            % Sum the resistance of each beam element
            sumMrBeams=sumMrBeams+Mp(elb,1);
        elseif nf(elb)==node
            % Sum the resistance of each beam element
            sumMrBeams=sumMrBeams+Mp(elb,2);
        end
    end
    MbNode(node,1)=sumMrBeams;
    
    % To detect which columns intersect with the node (if any)
    sumMrCols=0;
    for j=1:cols
        elc=elemcols(j);
        if ni(elc)==node || nf(elc)==node
            % Sum the resistance of each column element
            sumMrCols=sumMrCols+Mp(elc,1);
        end
    end
    McNode(node,1)=sumMrCols;
    EffMcb(node,1)=sumMrCols/sumMrBeams;
end
