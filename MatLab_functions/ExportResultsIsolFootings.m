function ExportResultsIsolFootings(directionData,bestDispositionFootings,...
    dimensionFootingCollection,nbarsFootingsCollection,typesRebarFooting,...
    coordBaseFooting,cols_sym_asym_isr)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsIsolFootings(directionData,bestDispositionFootings,...
%   dimensionFootingCollection,nbarsFootingsCollection,typesRebarFooting,...
%   coordBaseFooting,cols_sym_asym_isr)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of an isolated
% footing element into a .csv file on a prescribed folder route.
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
%                                       base, in format:
%
%                                               [x,y,z]
%
%                                       where z is the vertical axis
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
%         cols_sym_asym_isr:            is the parameter that indicates
%                                       which sort of reinforcement design
%                                       was performed: either symmetrical
%                                       rebar in columns, asymmetrical 
%                                       rebar in columns or the ISR
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-05
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
if isempty(directionData)==0
    % To export design results of isolated footings
    nombre_archivo='nrebar_footings_collection.csv';
    fileid_01=fopen([directionData,nombre_archivo],'w+t');

    nombre_archivo='dim_footings_collection.csv';
    fileid_02=fopen([directionData,nombre_archivo],'w+t');

    nombre_archivo='collection_type_bar_footings.csv';
    fileid_03=fopen([directionData,nombre_archivo],'w+t');

    nombre_archivo='collection_disposition_rebar_footings.csv';
    fileid_04=fopen([directionData,nombre_archivo],'w+t');

    % To wirte .csv file for the exportation of results
    nfoot=length(dimensionFootingCollection(:,1));
    if cols_sym_asym_isr~="ISR"
        fprintf(fileid_03,'%d\n',typesRebarFooting);

        for j=1:length(bestDispositionFootings(:,1))
            fprintf(fileid_04,'%.2f,%.2f\n',bestDispositionFootings(j,:));
        end          
        for i=1:nfoot
            fprintf(fileid_01,'%d,%d,%d,%d\n',nbarsFootingsCollection(i,:));
            fprintf(fileid_02,'%.2f,%.2f,%.2f,%.2f\n',...
                dimensionFootingCollection(i,:));
        end
    else
        for i=1:nfoot
            fprintf(fileid_02,'%.2f,%.2f,%.2f,%.2f\n',...
                dimensionFootingCollection(i,:));
        end
    end

    fclose(fileid_01);
    fclose(fileid_02);
    fclose(fileid_03);
    fclose(fileid_04);

    nombre_archivo='coord_base_footings.csv';

    fileid_05=fopen([directionData,nombre_archivo],'w+t');

    for i=1:nfoot
        fprintf(fileid_05,'%.2f,%.2f,%.2f\n',coordBaseFooting(i,:));
    end

    fclose(fileid_05);
end