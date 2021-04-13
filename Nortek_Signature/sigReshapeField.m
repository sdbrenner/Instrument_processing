function Data2D = sigReshapeField( Data, numSamples, fieldNames )
% SIGRESHAPREFIELD reshapes vector fields in data collected by Nortek
% Signature-Series ADCP's into 2D matrix forms
%
%	For data that are "well behaved", a simple 'reshape' is sufficient.
%	However, if bursts are cut off, this method fails.  In that case,
%	ensembles are delineated by indices where the ensemble count stops
%	increasing and re-organization is done in a loop
%
%   Data2D = sigReshapeField( Data, Config, fieldNames )



    ensembleCount = Data.Burst_EnsembleCount;

    % Loop through list of fields
    for m = 1:length(fieldNames)

        % Get data from specified field
        fieldNameM = fieldNames{m};
        fieldData = Data.(fieldNameM);

        % Check if field exists
        if ~isfield( Data, fieldNameM )
            warning('Incorrect field name, results replaced with NaN');
            Data2D.(fieldNameM) = NaN;
            continue;
        end


        % Reshape
        [~,numCol] = size( fieldData );
        try  % Try simple reshape:
            field2D = reshape( fieldData , numSamples,[],numCol );
        catch % If simple reshape fails, use loop:
            % Find indices when the ensemble count stops increasing:
            deltaCount = diff(double( ensembleCount ));
            negInd = find( deltaCount < 0 );
            negInd = [ 0; negInd; length( ensembleCount ) ];
            numEnsembles = length(negInd)-1; 

            % Pre-allocate "empty" (nan-filled) matrices
            field2D = NaN( numSamples, numEnsembles, numCol );

            % Loop through ensembles and fill matrices
            for n = 1:numEnsembles
                ensembleInd = (negInd(n)+1) : negInd(n+1) ;
                ensembleCountn = ensembleCount( ensembleInd );
                field2D( ensembleCountn, n, : ) =           fieldData( ensembleInd, : );
            end
        end

        % Save data in new struct
        Data2D.(fieldNameM) = field2D;
    end

end

%% EMBEDDED FUNCTIONS

% function field2D = reshape2D( fieldData, ensembleCount, numSamples)
% 
%     try  % Try simple reshape:
%         field2D = reshape( fieldData , numSamples,[] );
%     catch % If simple reshape fails, use loop:
%         % Find indices when the ensemble count stops increasing:
%         deltaCount = diff(double( ensembleCount ));
%         negInd = find( deltaCount < 0 );
%         negInd = [ 0; negInd; length( ensembleCount ) ];
%         numEnsembles = length(negInd)-1; 
% 
%         % Pre-allocate "empty" (nan-filled) matrices
%         field2D = NaN( numSamples, numEnsembles );
% 
%         % Loop through ensembles and fill matrices
%         for n = 1:numEnsembles
%             ensembleInd = (negInd(n)+1) : negInd(n+1) ;
%             ensembleCountn = ensembleCount( ensembleInd );
%             field2D( ensembleCountn, n ) =           fieldData( ensembleInd );
%         end
%     end
%     
% end
% 
% 
% 
% 
% function field2D = reshape3D( fieldData, ensembleCount, numSamples, numCol)
% 
%     try  % Try simple reshape:
%         field2D = reshape( fieldData , numSamples,[],numCol );
%     catch % If simple reshape fails, use loop:
%         % Find indices when the ensemble count stops increasing:
%         deltaCount = diff(double( ensembleCount ));
%         negInd = find( deltaCount < 0 );
%         negInd = [ 0; negInd; length( ensembleCount ) ];
%         numEnsembles = length(negInd)-1; 
% 
%         % Pre-allocate "empty" (nan-filled) matrices
%         field2D = NaN( numSamples, numEnsembles, numCol );
% 
%         % Loop through ensembles and fill matrices
%         for n = 1:numEnsembles
%             ensembleInd = (negInd(n)+1) : negInd(n+1) ;
%             ensembleCountn = ensembleCount( ensembleInd );
%             field2D( ensembleCountn, n, : ) =           fieldData( ensembleInd, : );
%         end
%     end
%     
% end




