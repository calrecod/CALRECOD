function ExportResultsColumn(directionData,dimColumnsCollection,...
    bestdisposicionRebar,nbarColumnsCollection,bestArrangement,...
    cols_sym_asym_isr,coordBaseCols)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsColumn(directionData,dimColumnsCollection,...
%  bestdisposicionRebar,nbarColumnsCollection,bestArrangement,coordBaseCols)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of a column element 
%          into a .csv file on a prescribed folder route.
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
%                                   cross-sections
%
%         bestArrangement:          is the list of the rebar diameters' 
%                                   indices used in the element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------

% To export design results 
nombre_archivo='coord_base_columns.csv';
fileid=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='nrebar_cols_collection.csv';
fileid_01=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='dim_cols_collection.csv';
fileid_02=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_type_bar_cols.csv';
fileid_03=fopen([directionData,nombre_archivo],'w+t');

nombre_archivo='collection_disposition_rebar_cols.csv';
fileid_04=fopen([directionData,nombre_archivo],'w+t');

% Writing the .csv file for the exportation of results
ncols=length(dimColumnsCollection(:,1));
if cols_sym_asym_isr~="ISR"
    for i=1:ncols
        fprintf(fileid_02,'%.2f,%.2f,%.2f,%.2f\n',dimColumnsCollection(i,:));
        fprintf(fileid_01,'%d\n',nbarColumnsCollection(i,:));
        fprintf(fileid,'%.2f,%.2f,%.2f\n',coordBaseCols(i,:));
    end
    fprintf(fileid_03,'%d\n',bestArrangement);

    for j=1:length(bestdisposicionRebar(:,1))
        fprintf(fileid_04,'%.2f,%.2f\n',bestdisposicionRebar(j,:));
    end
else
    for i=1:ncols
        fprintf(fileid_02,'%.2f,%.2f,%.2f,%.2f\n',dimColumnsCollection(i,:));
        fprintf(fileid,'%.2f,%.2f,%.2f\n',coordBaseCols(i,:));
    end
end
fclose(fileid);
fclose(fileid_01);
fclose(fileid_02);
fclose(fileid_03);
fclose(fileid_04);

