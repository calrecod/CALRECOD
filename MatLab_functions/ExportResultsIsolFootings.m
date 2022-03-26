function ExportResultsIsolFootings(directionData,bestDispositionFootings,...
    dimensionFootingCollection,nbarsFootingsCollection,typesRebarFooting,...
    coordBaseFooting)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsIsolFootings(directionData,bestDispositionFootings,...
%   dimensionFootingCollection,nbarsFootingsCollection,typesRebarFooting,...
%   coordBaseFooting)
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of an isolated
% footing element into a .txt file on a prescribed folder route.
%
%
% INPUT:  directionData:                is the folder disc location to save
%                                       the results
%
%         dimensionFootingCollection:   is the array containing the isolated
%                                       footing transversal cross-section
%                                       dimensions data
%
%         coordBaseFooting:             is the array containing the
%                                       coordinates of the isolated footing
%                                       base
%
%         bestDispositionFootings:      is the array containing the local
%                                       rebar coordinates of the isolated
%                                       footing transversal cross-sections
%
%         nbarsFootingsCollection:      is the total number of rebars on
%                                       the transversal cross-sections, both
%                                       in tension and compression
%
%         typesRebarFooting:            is the list of the rebar types used
%                                       in the element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

% To export design results __________________________________________
nombre_archivo='nrebar_footings_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_footings_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_footings.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_footings.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');

% To wirte .txt file for the exportation of results_______________
if cols_sym_asym_isr~="ISR"
    fprintf(fileid_03,'%d\n',typesRebarFooting);


    for j=1:length(bestDispositionFootings(:,1))
        fprintf(fileid_04,'%.2f %.2f\n',bestDispositionFootings(j,:));
    end          


    fprintf(fileid_01,'%d %d %d %d\n',nbarsFootingsCollection(i,:));


    fprintf(fileid_02,'%.2f %.2f %.2f %.2f\n',dimensionFootingCollection(i,:));

end

fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);
    
nombre_archivo='coord_desplante_footings.txt';

fileid_05=fopen([direction_data,nombre_archivo],'w+t');
for i=1:nzapatas
    fprintf(fileid_05,'%.2f %.2f %.2f\n',coordBaseFooting(i,:));
end
fclose(fileid_05);