function ExportResultsColumn(directionData,dimColumnsCollection,...
    bestdisposicionRebar,nbarColumnsCollection,bestArrangement,cols_sym_asym_isr,coordBaseCols)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsColumn(directionData,dimColumnsCollection,...
%  bestdisposicionRebar,nbarColumnsCollection,bestArrangement,coordBaseCols)
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of a column element 
%          into a .txt file on a prescribed folder route.
% 
% INPUT:  directionData:            is the folder disc location to save the
%                                   results
%
%         dimColumnsCollection:     is the array containing the cross-section
%                                   dimensions data of the column element

%         coordBaseCols:            is the array containing the coordinates 
%                                   of the column base cross-section's 
%                                   centroid
%
%         bestdisposicionRebar:     is the array containing the local rebar
%                                   coordinates of the column cross-sections
%
%         nbarColumnsCollection:    is the total number of rebars of column
%                                   cross-sections, both in tension and 
%                                   compression
%
%         bestArrangement:          is the list of the rebar types used in
%                                   the element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

% To export design results __________________________________________
nombre_archivo='nrebar_cols_collection.txt';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_cols_collection.txt';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_arerglo_bar_cols.txt';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_cols.txt';
fileid_04=fopen([directionData,nombre_archivo],'w+t');

% Writing the .txt file for the exportation of results______________
if cols_sym_asym_isr~="ISR"
    fprintf(fileid_02,'%.2f %.2f %.2f %.2f %.2f\n',dimColumnsCollection(1,:));

    fprintf(fileid_03,'%d\n',bestArrangement);

    for j=1:length(bestdisposicionRebar(:,1))
        fprintf(fileid_04,'%.2f %.2f\n',bestdisposicionRebar(j,:));
    end

    fprintf(fileid_01,'%d\n',nbarColumnsCollection(1));
end

fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

nombre_archivo='coord_desplante_columns.txt';

fileid=fopen([directionData,nombre_archivo],'w+t');
for j=1:cols
    fprintf(fileid,'%.2f %.2f %.2f\n',coordBaseCols(j,:));
end
fclose(fileid);
