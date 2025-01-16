% =========================================================================
% Calculate_precision =====================================================

addpath('src')

data_sets = {'small_data_set_out','large_data_set_out'};

for i = 1%:numel(data_sets)
   
    output_folder = fullfile('output',data_sets{i});
    mkSaveFolder(output_folder)   
    
    data_files = dir(fullfile('data',data_sets{i},'*trackPositions.csv'));
    
    precision = calculatePrecision(data_files);

    stats = calculateStatistics(precision,output_folder);

    % plotPrecision(precision,data_sets{i},output_folder);

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

    precision = [array2table(std_dev,'VariableNames',...
        {'x_std', 'y_std', 'z_std', 'overall'}), metaData];

end %----------------------------------------------------------------------


function metaData = parseMetadata(files) %----------------------------------
    
    if any(contains(files,'='))        

        mll2 = string(extractBetween(files, 'celltype=','_'));
        media = string(extractBetween(files, 'media=','_'));
        
        status = repmat("live",numel(files),1);
        fixed_idx = contains(files,'fixed','IgnoreCase', true);
        status(fixed_idx) = "fixed";        

        sample = strcat(string(extractBetween(files, 'date=','_')),...
            '_', string(extractBetween(files, 'well=','_')));
        replicate = extractBetween(files, 'Pos','_');
        
    else

        live_samples = ~contains(files,'fixed','IgnoreCase', true);
        
        files(live_samples) = strcat('LIVE_',files(live_samples));
        
        parsed_names = string(split(files,'_'));
        
        mll2 = upper(parsed_names(:,2));
        media = erase(parsed_names(:,3),'IF');      
        status = lower(parsed_names(:,1));
        sample = erase(parsed_names(:,5),'s');
        replicate = erase(parsed_names(:,6),'Pos');
        
    end

    metaData = table(mll2, media, status, sample, replicate);

end %----------------------------------------------------------------------


function std_dev = precisionCalculator(data_files) %-----------------------

    num_data_files = numel(data_files);
    
    std_dev = zeros(num_data_files,3);     
    
    for i = 1:num_data_files

       file = fullfile(data_files(i).folder,data_files(i).name);
       tracks = readtable(file,"VariableNamingRule","preserve");

        % Get center of mass for each track
        [mu, track_id, track_length] = grpstats(tracks{:,{'x','y','z'}}, ...
            tracks.("#track"), {'mean','gname', 'numel'});
    
        track_id = str2double(track_id);
        track_length = track_length(:, 1);
        
        % Move center of mass of all tracks to origin
        [~, idx] = ismember(tracks.("#track"),track_id);
        translated_tracks = tracks{:,{'x','y','z'}} - mu(idx,:);
    
        % remove all singleton trajectories
        idx = ismember(tracks.("#track"),track_id(track_length==1));
        translated_tracks(idx,:) = [];
    
        std_dev(i,:) = std(translated_tracks);
                          
    end

    overall = sqrt(sum(std_dev.^2,2));
    std_dev = [std_dev, overall];

end %----------------------------------------------------------------------


function stats = calculateStatistics(precision,output_folder) %------------

    group_names = ...
        strcat(precision.mll2,'_',precision.media,'_',precision.status);
    
    stats = struct('dim',{'x_std','y_std','z_std','overall'},...
        'anovaP',[],'anova',[],'anovaMulticomp',[],'anovaGnames',[],...
        'kwP',[],'kw',[],'kwMulticomp',[],'kwGnames',[]);

    if contains(output_folder,'small')
        CritValType = "lsd";
    else
        CritValType = "hsd";
    end

    for i = 1:numel(stats)

        % N-Way ANOVA (parametric)
        [stats(i).anovaP,~,stats(i).anova] = ...
            anovan(precision.(stats(i).dim),{group_names},"display","off");
        [stats(i).anovaMulticomp,~,~,stats(i).anovaGnames] = ...
            multcompare(stats(i).anova,"CriticalValueType",CritValType,"display","off");        

        % Kruskal-Wallis Test (nonparametric)
        [stats(i).kwP,~,stats(i).kw] = kruskalwallis(...
            precision.(stats(i).dim),group_names,"off");
        [stats(i).kwMulticomp,~,~,stats(i).kwGnames] = ...
            multcompare(stats(i).kw,"CriticalValueType",CritValType,"display","off");                   

    end

    saveStats(stats);

    save(fullfile(output_folder,"stats"),"stats")

end %----------------------------------------------------------------------


function saveStats(stats) %------------------------------------------------
    
    
    
end %----------------------------------------------------------------------


function plotPrecision(precision,output_folder) %--------------------------

    group_names = ...
        strcat(precision.mll2,'_',precision.media,'_',precision.status);

    figure(1)
    clf

    if numel(group_names) == 5

    else
    
    end
    boxplot(precision.x_std,...
        strcat(precision.mll2,'_',precision.media,'_',precision.status));

    ylim([0,max(precision.x_std)])


end %----------------------------------------------------------------------

% FUNCTIONS ===============================================================
% =========================================================================