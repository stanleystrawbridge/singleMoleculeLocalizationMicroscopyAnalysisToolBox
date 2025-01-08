% =========================================================================
% Calculate_precision =====================================================

addpath('src')


%% Small Data Set ---------------------------------------------------------

mkSaveFolder(save_folder)   
data_files = ...
    dir(fullfile('data','small_data_set_out','*trackPositions.csv'));

precision = calculatePrecision(data_files);




%% Large Data Set ---------------------------------------------------------
data_files = ...
    dir(fullfile('data','large_data_set_out','*trackPositions.csv'));

calculatePrecision(data_files);

% Calculate_precision =====================================================
% =========================================================================


% =========================================================================
% FUNCTIONS ===============================================================

function mkSaveFolder(save_folder)    
    
    if ~exist(save_folder,'dir')
        mkdir(save_folder)
    end

end %----------------------------------------------------------------------


function precision = calculatePrecision(data_files) %----------------------
    
    num_data_files = numel(data_files);
    
    prceision = struct('std',zeros(num_data_files,1),...
        'cell',...
        'media',...
        'status',
        'date',...        
        'well'
        'replicate',...   )

    for i = 1:num_data_files

           file = fullfile(data_files(i).folder,data_files(i).name);
           tracks = readtable(file,"VariableNamingRule","preserve");

           stdDev = precisionCalculator(tracks);

           metaData = parseMetadata(data_files(i).name);
                       
    end

end %----------------------------------------------------------------------


function stdDev = precisionCalculator(tracks) %-------------------------

    % Get center of mass for each track
    [mu, track_id, track_length] = ...
        grpstats(tracks{:,{'x','y','z'}}, tracks.("#track"), ...
        {'mean','gname', 'numel'});

    track_id = str2double(track_id);
    track_length = track_length(:, 1);
    
    % Move center of mass of all tracks to origin
    [~, idx] = ismember(tracks.("#track"),track_id);

    translated_tracks = tracks{:,{'x','y','z'}} - mu(idx,:);

    % remove all singleton trajectories
    idx = ismember(tracks.("#track"),track_id(track_length==1));
    translated_tracks(idx,:) = [];

    stdDev = std(translated_tracks);
    
end %----------------------------------------------------------------------


function parseMetadata(file) %---------------------------------------------

    if contains(file,'=')
        
    else
        
    end

end %----------------------------------------------------------------------


function plotPrecision(precision) %----------------------------------------

    

end %----------------------------------------------------------------------

% FUNCTIONS ===============================================================
% =========================================================================
