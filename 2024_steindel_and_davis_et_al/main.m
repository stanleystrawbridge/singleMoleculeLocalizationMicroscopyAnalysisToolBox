% =========================================================================
% Calculate_precision =====================================================

addpath('src')

data_sets = {'small_data_set_out','large_data_set_out'};

for i = 1:numel(data_sets)
   
    output_folder = fullfile('output',data_sets{i});
    mkSaveFolder(output_folder)   
    
    data_files = dir(fullfile('data',data_sets{i},'*trackPositions.csv'));
    
    precision = calculatePrecision(data_files);
    
    % save precision

    % plotPrecision(precision);

end

% Calculate_precision =====================================================
% =========================================================================


% =========================================================================
% FUNCTIONS ===============================================================

function mkSaveFolder(output_folder)    
    
    if ~exist(output_folder,'dir')
        mkdir(output_folder)
    end

end %----------------------------------------------------------------------


function precision = calculatePrecision(data_files) %----------------------
    
    metaData = parseMetadata({data_files.name}');

    std_dev = precisionCalculator(data_files);

    precision = join(std_dev, metaData);

end %----------------------------------------------------------------------


function metaData = parseMetadata(files) %----------------------------------

    
    if any(contains(files,'='))
    
        mll2 = strings(num_data_files,1);
        media = strings(num_data_files,1);
        status = strings(num_data_files,1);
        sample = strings(num_data_files,1);
        replicate = strings(num_data_files,1);
        
    else

        fixed_samples = contains(files,'fixed','IgnoreCase', true);

        mll2 = strings(num_data_files,1);
        media = strings(num_data_files,1);
        status = strings(num_data_files,1);
        sample = strings(num_data_files,1);
        replicate = strings(num_data_files,1);
        
    end

    metaData = table(mll2

end %----------------------------------------------------------------------


function std_dev = precisionCalculator(data_files) %-----------------------

    num_data_files = numel(data_files);
    
    std_dev = zeros(num_data_files,1);     
    
    for i = 1:num_data_files

       file = fullfile(data_files(i).folder,data_files(i).name);
       tracks = readtable(file,"VariableNamingRule","preserve");

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
    
        std_dev(i) = std(translated_tracks);
                          
    end

end %----------------------------------------------------------------------


function plotPrecision(precision) %----------------------------------------

    

end %----------------------------------------------------------------------

% FUNCTIONS ===============================================================
% =========================================================================
