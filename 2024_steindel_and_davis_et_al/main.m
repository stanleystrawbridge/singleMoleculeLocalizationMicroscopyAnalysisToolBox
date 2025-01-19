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
    
    % plotPrecision(precision,stats,output_folder);

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
    unique_groups = unique(group_names);

    stats.precision_dim = string(precision.Properties.VariableNames(1:4));
    N_precision_types = numel(stats.precision_dim);    

    if contains(output_folder,'small')        
        
        % compare all live groups to fixed
        fixed_sample_name_idx = contains(unique_groups,'fixed');
        fixed_group_name = unique_groups(fixed_sample_name_idx);

        fixed_sample_idx = strcmp(group_names,fixed_group_name);
        fixed_sample = precision(fixed_sample_idx,:);

        live_group_names = unique_groups(~fixed_sample_name_idx);            
        N_live_groups = numel(live_group_names);

        stats.group_name = live_group_names;
        stats.p_value = zeros(N_live_groups, N_precision_types);               
        
        stats.live_mu = stats.p_value;
        stats.live_sigma = stats.p_value;
        stats.live_n =  stats.p_value(:,1);

        stats.fixed_mu = mean(fixed_sample{:,1:4});
        stats.fixed_sigma = std(fixed_sample{:,1:4});
        stats.fixed_n = size(fixed_sample{:,1:4},1);

        for i = 1:N_live_groups

            live_sample_idx = strcmp(group_names,live_group_names(i));
            live_sample = precision(live_sample_idx,:);
        
            stats.live_mu(i,:) = mean(live_sample{:,1:4});
            stats.live_sigma(i,:) = std(live_sample{:,1:4});
            stats.live_n(i) =  size(live_sample{:,1:4},1);

            for j = 1:N_precision_types
                
                stats.p_value(i,j) = ranksum(...
                    fixed_sample.(stats.precision_dim(j)),...
                    live_sample.(stats.precision_dim(j)),...
                    "tail","left");

            end
    
        end
        
        % Bonferoni correction
        stats.p_value = stats.p_value * N_live_groups;

    else

        % compare live and fixed camples on a genotype/media combo basis        
        stats.group_name = erase(unique_groups,'_fixed');
        stats.group_name = erase(stats.group_name,'_live');
        stats.group_name = unique(stats.group_name);
        N_groups = numel(fixed_group_names);

        stats.p_value = zeros(N_groups, N_precision_types);

        stats.live_mu = stats.p_value;
        stats.live_sigma = stats.p_value;
        stats.live_n =  stats.p_value(:,1);

        stats.fixed_mu = stats.p_value;
        stats.fixed_sigma = stats.p_value;
        stats.fixed_n =  stats.p_value(:,1);

        for i = 1:N_groups            
            
            fixed_sample_idx = strcmp(group_names,...
                strcat(stats.group_name(i), '_fixed'));
            fixed_sample = precision(fixed_sample_idx,:);

            live_sample_idx = strcmp(group_names,...
                strcat(stats.group_name(i),'_live'));
            live_sample = precision(live_sample_idx,:);

            stats.live_mu(i,:) = mean(live_sample{:,1:4});
            stats.live_sigma(i,:) = std(live_sample{:,1:4});
            stats.live_n(i) = size(live_sample{:,1:4},1);
    
            stats.fixed_mu(i,:) = mean(fixed_sample{:,1:4});
            stats.fixed_sigma(i,:) = std(fixed_sample{:,1:4});
            stats.fixed_n(i) =  size(fixed_sample{:,1:4},1);

            for j = 1:N_precision_types
            
                stats.p_value(i,j) = ranksum(...
                    fixed_sample.(stats.precision_dim(j)),...
                    live_sample.(stats.precision_dim(j)),...
                    "tail","left");
            
            end

        end

    end    

    % saveStats(stats,output_folder)

end %----------------------------------------------------------------------


function saveStats(stats,output_folder) %----------------------------------

    % Write data to m file_________________________________________________
    obj_save_name = fullfile(output_folder,"obj_stats");
    save(obj_save_name,"stats")
    
    % Write data to text file______________________________________________

    txt_save_name = fullfile(output_folder,"stats.txt");
    fileID = fopen( txt_save_name,'w'); 

    fprintf(fileID,'%-20s\n',...
        string(extractBetween(output_folder, '\','_out')));
    fprintf(fileID,'\n\n');

    if contains(output_folder,'small')

        

    else



    end
    
    fclose(fileID);

    
    
       
    
   

    fprintf(fileID,'GROUPS\n');    
    for i = 1:numel(stats(1).kwGnames)
        fprintf(fileID,'%-2s %-10s\n', num2str(i), stats(1).kwGnames{i});
    end
    fprintf(fileID,'\n\n');

    for i = 1:numel(stats)    
        fprintf(fileID,'precision %-10s\n\n',stats(i).dim);    
        fprintf(fileID,'%-4s %-4s %-7s %-7s \n','g1','g2','ANOVA','K-W');
        fprintf(fileID,'%-4u %-4u %-7.4f %-7.4f \n',...
            [stats(i).anovaMulticomp(:,[1 2 6]),...
            stats(i).kwMulticomp(:,6)]');
        fprintf(fileID,'\n\n');
    end


end %----------------------------------------------------------------------


function plotPrecision(precision,stats,output_folder) %--------------------

    if contains(output_folder,'small')

        

    else



    end

    % group_names = ...
    %     strcat(precision.mll2,'_',precision.media,'_',precision.status);
    % 
    % figure(1)
    % clf
    % 
    % if numel(group_names) == 5
    % 
    % else
    % 
    % end
    % boxplot(precision.x_std,...
    %     strcat(precision.mll2,'_',precision.media,'_',precision.status));
    % 
    % ylim([0,max(precision.x_std)])

end %----------------------------------------------------------------------

% FUNCTIONS ===============================================================
% =========================================================================




% ARCHIVE

% CritValType = "lsd";
%
% group_names = ...
%     strcat(precision.mll2,'_',precision.media,'_',precision.status);
% unique_groups = unique(group_names);
% 
% stats = struct('dim',{'x_std','y_std','z_std','overall'},...
%     'anovaP',[],'anova',[],'anovaMulticomp',[],'anovaGnames',[],...
%     'kwP',[],'kw',[],'kwMulticomp',[],'kwGnames',[]); 
% for i = 1:numel(stats)
% 
%     % N-Way ANOVA (parametric)
%     [stats(i).anovaP,~,stats(i).anova] = ...
%         anovan(precision.(stats(i).dim),{group_names},"display","off");
%     [stats(i).anovaMulticomp,~,~,stats(i).anovaGnames] = ...
%         multcompare(stats(i).anova,"CriticalValueType",...
%         CritValType,"display","off");        
% 
%     % Kruskal-Wallis Test (nonparametric)
%     [stats(i).kwP,~,stats(i).kw] = kruskalwallis(...
%         precision.(stats(i).dim),group_names,"off");
%     [stats(i).kwMulticomp,~,~,stats(i).kwGnames] = ...
%         multcompare(stats(i).kw,"CriticalValueType",...
%         CritValType,"display","off");                   
% 
% end        
% 
% saveStats(stats,output_folder);