function [VcrxNode,VuxNode,VxSF,hnodeMin]=shearNodes2DFramesSI(ni,nf,...
 elembeams,elemcols,fcNodes,fy,MrBeamsCollec,steelAreaBeams3sec,...
 dimensions,nnodes,coordBaseCols,L)
%-------------------------------------------------------------------------
% SYNTAX : 
% [VcrxNode,VuxNode,VxSF,hnodeMin]=shearNodes2DFramesSI(ni,nf,elembeams,...
% elemcols,Mp,fcNodes,fy,MrBeamsCollec,steelAreaBeams3sec,dimensions,...
% nnodes,coordBaseCols,L)
%-------------------------------------------------------------------------
%    PURPOSE
%     To revise the shear resistance of a Reinforced Concrete 2D Frame's 
%     nodes in the International System of Units, (kg,cm).
% 
%    OUTPUT:  hnodeMin:                 is the minimum node height
%                                       dimension to withstand the shear
%                                       demand Vux. Size: [nnodes, 1]
%
%             VxSF:                     is the Shear Efficiency of each 
%                                       node in the global X direction, 
%                                       computed as Vu/Vcr.
%
%             VcrxNode:                 resistant shear of the node
%
%             VuxNode:                  shear demand for the node
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
%             fcNodes:                  is the f'c of the nodes
%
%             fy:                       is the yield stress of the
%                                       reinforcing steel
%
%             MrBeamsCollec:            is the list of resistant bending
%                                       moments for three cross-sections 
%                                       of each beam element.
%                                           Size: n-beams x 3
%
%             steelAreaBeams3sec:       is the quantity of reinforcing area
%                                       of the three cross-sections of each
%                                       beam. The reinforcing steel area in
%                                       tension and compression is listed
%                                       for each cross-section, therefore:
%                                           Size: n-beams x 6
%
%                                           Format: 
%                                [area-tension-left, area-comp-left,...
%                                 area-tension-center,area-comp-center,...
%                                 area-tension-right, area-comp-right]
%
%              dimensions:              cross-section dimensions of each
%                                       element. Size: n-elems x 2
%
%                                                Format: [width, height]
%
%              nnodes:                  number of nodes of the 2D Frame
%
%              coordBaseCols:           list of coordinates of the base of
%                                       each column element in the 3D space
%                                       of reference, in which the height
%                                       dimension corresponds to the Z
%                                       axis. Size: n-cols x 3
%
%                                             Format: [x, y, z]
%
%               L:                      length of each element. 
%                                           Size: n-elems x 1
%
%-------------------------------------------------------------------------
% LAST MODIFIED: L.Verduzco    2023-07-03
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%-------------------------------------------------------------------------
beams=length(elembeams);
cols=length(elemcols);
steelAreaBeamsNode=[steelAreaBeams3sec(:,1:2),steelAreaBeams3sec(:,5:6)];
MrBeamsNode=[MrBeamsCollec(:,1), MrBeamsCollec(:,3)];

VcrxNode=zeros(nnodes,1);
VuxNode=zeros(nnodes,1);
VxSF=zeros(nnodes,1);
hnodeMin=zeros(nnodes,1);

% To detect the initial and final node of each beam element
for i=1:beams
    NiNfBeam(i,1)=ni(elembeams(i));
    NiNfBeam(i,2)=nf(elembeams(i));
end

% REVISION AND DESIGN OF NODES AGAINST SHEAR
for node=1:nnodes
    % To detect which beam element intersect with the node (if any)
    beamsNode=[];
    nb=0;
    xbeams=[];
    for j=1:beams
        
        elb=elembeams(j);
        if NiNfBeam(j,1)==node
            nb=nb+1;
            xbeams=[xbeams;
                    nb];
            beamsNode=[beamsNode;
                        j,1];
        elseif NiNfBeam(j,2)==node
            nb=nb+1;
            xbeams=[xbeams;
                    nb];
                
            beamsNode=[beamsNode;
                        j,2];
        end
    end
    % To detect which columns intersect with the node (if any)
    colNode=[];
    for j=1:cols
        elc=elemcols(j);
        if ni(elemcols(j))==node || nf(elemcols(j))==node
            colNode=[colNode;
                      j];
        end
    end
    if isempty(beamsNode)==0
        nbeamNode=length(beamsNode(:,1));
    else
        nbeamNode=[0];
    end
    ncolNode=length(colNode);
    if nbeamNode==1 && ncolNode==2
        %% EXTREME NODE - TYPE 2
        % ---------------------------------------------------------------
        % Shear in the global X direction
        % ---------------------------------------------------------------
        
        sb1=beamsNode(xbeams(1),2);

        % Calculate the normal force that the node withstands by the
        % beam element
        Fvx=fy*(steelAreaBeamsNode(beamsNode(xbeams(1),1),sb1*2-1))*1.25;

        Vcolx=abs(MrBeamsNode(beamsNode(xbeams(1),1),sb1))/...
            (0.5*L(elemcols(colNode(1)))+...
             0.5*L(elemcols(colNode(2))));

        % Node's effective width
        bv=dimensions(elembeams(beamsNode(xbeams(1),1)),1);

        if coordBaseCols(colNode(1),3)>coordBaseCols(colNode(2),3)

            % if the column 1 is superior
            bc=dimensions(elemcols(colNode(1)),1);
            hc=dimensions(elemcols(colNode(1)),2);
        elseif coordBaseCols(colNode(1),3)<coordBaseCols(colNode(2),3)

            % if the column 2 is inferior
            bc=dimensions(elemcols(colNode(2)),1);
            hc=dimensions(elemcols(colNode(2)),2);
        end
        be1=(bv+bc)*0.5;
        be2=bv+hc;
        be3=bc;

        be=min([be1,be2,be3]);

        Vux=Fvx-Vcolx;

        Vcrx=4.5*0.75*sqrt(fcNodes)*be*hc;
        VcrxNode(node,1)=Vcrx;
        VuxNode(node,1)=Vux;
        VxSF(node,1)=Vux/Vcrx;
        
        hnodeMin(node,1)=Vux/(4.5*0.75*sqrt(fcNodes)*be);
        
    elseif nbeamNode==2 && ncolNode==2
        %% INTERIOR NODE - TYPE 2
        % -------------------------------------------------------------
        % Shear in the global X direction
        % -------------------------------------------------------------
        sb1=beamsNode(xbeams(1),2);
        sb2=beamsNode(xbeams(2),2);

        % Calculate the normal force that the node withstands by the
        % beam element
        Fvx=fy*(steelAreaBeamsNode(beamsNode(xbeams(1),1),sb1*2-1)+...
            steelAreaBeamsNode(beamsNode(xbeams(2),1),sb2*2))*1.25;

        % Determine the shear contribution by the column
        if MrBeamsNode(beamsNode(xbeams(1),1),sb1)>...
                MrBeamsNode(beamsNode(xbeams(2),1),sb2)

            Vcolx=abs(MrBeamsNode(beamsNode(xbeams(1),1),sb1))/...
                (0.5*L(elemcols(colNode(1)))+...
                 0.5*L(elemcols(colNode(2))));
        elseif MrBeamsNode(beamsNode(xbeams(1),1),sb1)<=...
                MrBeamsNode(beamsNode(xbeams(2),1),sb2)

            Vcolx=abs(MrBeamsNode(beamsNode(xbeams(2),1),sb2))/...
                (0.5*L(elemcols(colNode(1)))+...
                 0.5*L(elemcols(colNode(2))));
        end
        % Node's effective width
        bv=0.5*(dimensions(elembeams(beamsNode(xbeams(1),1)),1)+...
                dimensions(elembeams(beamsNode(xbeams(2),1)),1));
        
        if coordBaseCols(colNode(1),3)>coordBaseCols(colNode(2),3)

            % if the column 1 is superior
            bc=dimensions(elemcols(colNode(1)),1);
            hc=dimensions(elemcols(colNode(1)),2);
        elseif coordBaseCols(colNode(1),3)<coordBaseCols(colNode(2),3)

            % if the column 2 is inferior
            bc=dimensions(elemcols(colNode(2)),1);
            hc=dimensions(elemcols(colNode(2)),2);
        end
        be1=(bv+bc)*0.5;
        be2=bv+hc;
        be3=bc;

        be=min([be1,be2,be3]);

        Vux=Fvx-Vcolx;

        Vcrx=5.5*0.75*sqrt(fcNodes)*be*hc;
        VcrxNode(node,1)=Vcrx;
        VuxNode(node,1)=Vux;
        VxSF(node,1)=Vux/Vcrx;
        hnodeMin(node,1)=Vux/(5.5*0.75*sqrt(fcNodes)*be);
        
    elseif nbeamNode==2 && ncolNode==1
        %% INTERIOR ROOF NODE - TYPE 2
        % -------------------------------------------------------------
        % Shear revision in the global X direction
        % -------------------------------------------------------------
        sb1=beamsNode(xbeams(1),2);
        sb2=beamsNode(xbeams(2),2);

        % Calculate the normal force that the node withstands by the
        % beam element
        Fvx=fy*(steelAreaBeamsNode(beamsNode(xbeams(1),1),sb1*2-1)+...
            steelAreaBeamsNode(beamsNode(xbeams(2),1),sb2*2))*1.25;

        % Determine the shear action provided by the column
        Vcolx=0;

        % Node's effective width
        bv=0.5*(dimensions(elembeams(beamsNode(xbeams(1),1)),1)+...
                dimensions(elembeams(beamsNode(xbeams(2),1)),1));
            
        bc=dimensions(elemcols(colNode(1)),1);
        hc=dimensions(elemcols(colNode(1)),2);

        be1=(bv+bc)*0.5;
        be2=bv+hc;
        be3=bc;

        be=min([be1,be2,be3]);

        Vux=Fvx-Vcolx;

        Vcrx=5.5*0.75*sqrt(fcNodes)*be*hc;
        
        VcrxNode(node,1)=Vcrx;
        VuxNode(node,1)=Vux;
        VxSF(node,1)=Vux/Vcrx;
        hnodeMin(node,1)=Vux/(5.5*0.75*sqrt(fcNodes)*be);
        
    elseif nbeamNode==1 && ncolNode==1
        %% EXTERIOR ROOF NODE - TYPE 2 (CORNER)
        % -------------------------------------------------------------
        % Shear revision in the global X direction
        % -------------------------------------------------------------
        sb1=beamsNode(xbeams(1),2);
        % Calculate the normal force that the node withstands by the
        % beam element
        Fvx=fy*(steelAreaBeamsNode(beamsNode(xbeams(1),1),sb1*2-1))*1.25;
        Vcolx=0;
        
        % Node's effective width
        bv=dimensions(elembeams(beamsNode(xbeams(1),1)),1);
            
        bc=dimensions(elemcols(colNode(1)),1);
        hc=dimensions(elemcols(colNode(1)),2);

        be1=(bv+bc)*0.5;
        be2=bv+hc;
        be3=bc;

        be=min([be1,be2,be3]);

        Vux=Fvx-Vcolx;

        Vcrx=3.5*0.75*sqrt(fcNodes)*be*hc;
        VcrxNode(node,1)=Vcrx;
        VuxNode(node,1)=Vux;
        VxSF(node,1)=Vux/Vcrx;
        hnodeMin(node,1)=Vux/(3.5*0.75*sqrt(fcNodes)*be);
        
    end
end