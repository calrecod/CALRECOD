function ExportResultsBeam(directionData,dim_beams_collection,coordEndBeams,...
    disposition_rebar_beams3sec,nbarbeamsCollection,arrangemetbarbeams)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsBeam(directionData,dim_beams_collection,coordEndBeams,...
%   disposition_rebar_beams3sec,nbarbeamsCollection,arrangemetbarbeams)
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of a beam element
%          into a .txt file on a prescribed folder route.
% 
% INPUT:  directionData:            is the folder disc location to save the
%                                   results
%
%         dim_beams_collection:     is the array containing the cross-section
%                                   dimensions data of the beam element
%
%         coordEndBeams:            is the array containing the coordinates 
%                                   of the initial end's cross-section 
%                                   centroid of the beam
%
%         disposition_rebar_beams3sec: is the array containing the local 
%                                      rebar coordinates of each of the 
%                                      three critical design cross-sections
%                                      of the beam
%
%         nbarbeamsCollection:      is the total number of rebars of each 
%                                   of the three design cross-sections, both 
%                                   in tension and compression. 
%                                   Size = [1,6] in format:
%
%               [nbarsLeft_{ten},nbarsLeft_{comp},nbarsCenter_{ten},...
%               nbarsCenter_{comp},nbarsRight_{ten},nbarsRight_{comp}]
%
%         arrangemetbarbeams:       is the list of all the rebar types used 
%                                   in the element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

% To export design results ___________________________________________
nombre_archivo='nrebar_beams_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_beams_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_beams.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_beams3sec.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');

% To write .txt files for the exportation of results _____________
if cols_sym_asym_isr~="ISR"
    fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',...
                    dim_beams_collection(i,:));


    for j=1:length(disposition_rebar_beams3sec(:,1))
        fprintf(fileid_04,'%.2f %.2f\n',disposition_rebar_beams3sec(j,:));
    end

    fprintf(fileid_01,'%d %d %d %d %d %d\n',nbarbeamsCollection(i,:));

    fprintf(fileid_03,'%d\n',arrangemetbarbeams(:,:));
else
    fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',...
                    dim_beams_collection(i,:));

    fprintf(fileid_04,'%.2f %.2f\n',[0,0]);


    fprintf(fileid_01,'%d %d %d %d %d %d\n',[0,0,0,0,0,0]);

    fprintf(fileid_03,'%d\n',[0]);
end
fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

nombre_archivo='coord_comienzo_beams.txt';
fileid=fopen([direction_data,nombre_archivo],'w+t');
for i=1:beams
    fprintf(fileid,'%.2f %.2f %.2f\n',coordEndBeams(i,:));
end
fclose(fileid);
