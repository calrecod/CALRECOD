function [dimensions]=AdjustDimElemFrames(nbars,type_elem,...
    nnodes,dimensions,coord_base_cols,ni,nf)

%------------------------------------------------------------------------
% Syntax:
% [dimensions]=AdjustDimElemFrames(nbars,type_elem,...
%   nnodes,dimensions,coord_base_cols,ni,nf)
%
%------------------------------------------------------------------------
% PURPOSE: To adjust cross-section dimensions of elements accordingly 
% to make them uniform.
% 
% OUTPUT: dimensions:               are the new modified cross-section
%                                   dimensions
%
% INPUT:  nbars:                    number of elements composing the 
%                                   structural frame
%
%         type_elem:                is the identification vector containing
%                                   the label for each elements as ''Beam''
%                                   or ''Col'' for each of the elements.
%                                   Size [nbars,2] in format [nelem,label]
%
%         nnodes:                   total number of nodes of the structural
%                                   frame to be reviewed
%
%         dimensions:               is the array containing the cross-
%                                   section dimensions data for each of the
%                                   elements. Size = [nbars,2] in format 
%                                   [be,he]
%
%         coord_base_cols:          is the array containing the coordinates 
%                                   of the lower ends' centroid of each 
%                                   column. Size = [ncolumns,2] in format 
%                                   [xc,yc]
%
%         ni,nf:                    vectors containing the ID of initial 
%                                   nodes and final nodes of each element, 
%                                   respectively, in the order 
%                                   pre-established
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%________________________________________________________________________
%%% UPDATE FOR THE DIMENSIONS OF ELEMENTS - UNIFORMITY OF WIDTH____
%________________________________________________________________________
elem_cols=[];
elem_beams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elem_cols=[elem_cols,j];
    end
end


ni_nf_beam=zeros(beams,2);
for i=1:beams
    ni_nf_beam(i,1)=ni(elem_beams(i));
    ni_nf_beam(i,2)=nf(elem_beams(i));
end

for node=1:nnodes
    % To detect if there are beams converging to that node in question__
    beams_node=[];
    for j=1:beams
        if ni_nf_beam(j,1)==node
            beams_node=[beams_node;
                            j,1];
        elseif ni_nf_beam(j,2)==node
            beams_node=[beams_node;
                            j,1];
        end
    end
    
    % To detect if there are columns converging to that node in question__
    col_node=[];
    for j=1:cols
        if ni(elem_cols(j))==node || nf(elem_cols(j))==node
            col_node=[col_node;
                      j];

        end
    end
    
    
    if isempty(beams_node)==1
            % Support node of column-footing___________________________
    elseif isempty(beams_node)==0
        if length(beams_node(:,1))==1 && length(col_node)==2
            % Inter-story Extreme node
            bv=dimensions(elem_beams(beams_node(1,1)),1);
            dimensions(elem_cols(col_node(2)),1);
            dimensions(elem_cols(col_node(1)),1);
            
            if bv>dimensions(elem_cols(col_node(2)),1) && bv>dimensions(elem_cols(col_node(1)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
                dimensions(elem_cols(col_node(2)),1)=bv;
                
            elseif bv>dimensions(elem_cols(col_node(1)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
                
            elseif bv>dimensions(elem_cols(col_node(2)),1)
                dimensions(elem_cols(col_node(2)),1)=bv;
                
            end
            dimensions(elem_cols(col_node(1)),1);
            dimensions(elem_cols(col_node(2)),1);
            
            if coord_base_cols(col_node(1),2)>coord_base_cols(col_node(2),2)
                % IF COL 1 IS ABOVE __________________________________
                if dimensions(elem_cols(col_node(2)),1)<dimensions(elem_cols(col_node(1)),1)
                    dimensions(elem_cols(col_node(2)),1)=dimensions(elem_cols(col_node(1)),1);
                end
                
            elseif coord_base_cols(col_node(1),2)<coord_base_cols(col_node(2),2)
                % IF COL 2 IS ABOVE __________________________________
                if dimensions(elem_cols(col_node(1)),1)<dimensions(elem_cols(col_node(2)),1)
                    dimensions(elem_cols(col_node(1)),1)=dimensions(elem_cols(col_node(2)),1);
                end
                
                
            end
            
        elseif length(beams_node(:,1))==2 && length(col_node)==2
            % Inter-story interior node__________________________________
            bv_01=dimensions(elem_beams(beams_node(1,1)),1);
            bv_02=dimensions(elem_beams(beams_node(2,1)),1);

            bv=max([bv_01, bv_02]);
            dimensions(col_node(2),1);
            dimensions(col_node(1),1);
            
            if bv>dimensions(elem_cols(col_node(1)),1) && bv>dimensions(elem_cols(col_node(2)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
                dimensions(elem_cols(col_node(2)),1)=bv;
            elseif bv>dimensions(elem_cols(col_node(1)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
            elseif bv>dimensions(elem_cols(col_node(2)),1)
                dimensions(elem_cols(col_node(2)),1)=bv;
            end
            dimensions(elem_cols(col_node(2)),1);
            dimensions(elem_cols(col_node(1)),1);
            
            if coord_base_cols(col_node(1),2)>coord_base_cols(col_node(2),2)
                % IF COL 1 IS ABOVE ________________________________
                if dimensions(elem_cols(col_node(2)),1)<dimensions(elem_cols(col_node(1)),1)
                    dimensions(elem_cols(col_node(2)),1)=dimensions(elem_cols(col_node(1)),1);
               
                end
                
                
            elseif coord_base_cols(col_node(1),2)<coord_base_cols(col_node(2),2)
                % IF COL 2 IS ABOVE ________________________________
                if dimensions(elem_cols(col_node(1)),1)<dimensions(elem_cols(col_node(2)),1)
                    dimensions(elem_cols(col_node(1)),1)=dimensions(elem_cols(col_node(2)),1);
                end
                
                
            end
            
            
        elseif length(beams_node(:,1))==2 && length(col_node)==1
            %fprintf("\nNodo interior azotea o planta baja\n");
            
            bv_01=dimensions(elem_beams(beams_node(1,1)),1);
            bv_02=dimensions(elem_beams(beams_node(2,1)),1);
            bv=max([bv_01, bv_02]);
            dimensions(elem_cols(col_node(1)),1);
            if bv>dimensions(elem_cols(col_node(1)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
            end
            dimensions(elem_cols(col_node(1)),1);
            
        elseif length(beams_node(:,1))==1 && length(col_node)==1
            % Rooftop corner node_________________________________
            bv=dimensions(elem_beams(beams_node(1,1)),1);
            dimensions(elem_cols(col_node(1)),1);
            if bv>dimensions(elem_cols(col_node(1)),1)
                dimensions(elem_cols(col_node(1)),1)=bv;
            end
            dimensions(elem_cols(col_node(1)),1);
        
        end
        
    end
end

% TO UNIFORM COLUMN ALONG AN AXIS ______________________________________
% Note: The lower columns will have cross-section dimensions equal or
% bigger than the upper ones ...........................................
for node=1:nnodes

    % To detect if there are columns intersecting the node in question___
    col_node=[];
    for j=1:cols
        if ni(elem_cols(j))==node || nf(elem_cols(j))==node
            col_node=[col_node;
                      j];

        end
    end

    % To verify that the upper columns have lower or equal height
    % dimensions than the lower ones ____________________________________
    if length(col_node)==2
        if coord_base_cols(col_node(1),2)>coord_base_cols(col_node(2),2)
            % IF COL 1 IS ABOVE
            if dimensions(elem_cols(col_node(2)),2)<dimensions(elem_cols(col_node(1)),2)
                dimensions(elem_cols(col_node(2)),2)=dimensions(elem_cols(col_node(1)),2);

            end

        elseif coord_base_cols(col_node(2),2)>coord_base_cols(col_node(1),2)
            % IF COL 2 IS ABOVE ______________________________________
            if dimensions(elem_cols(col_node(1)),2)<dimensions(elem_cols(col_node(2)),2)
                dimensions(elem_cols(col_node(1)),2)=dimensions(elem_cols(col_node(2)),2);

            end
        end
    end

end  

