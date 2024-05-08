% Cropping for Double Helix Point Spead Function SMLM Data
% 1.) Place all coordinate csv's into a single data_path folder
% 2.) provide the full data path in user inputs (E.g. C:\my_data_folder)
% 3.) Set wether data has tracked or not (true or false)
% 4.) Set the limits for cropping in nm, usally +/- 2000nm. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path = 'D:\test\raw_data';
tracked = false; % true/flase
z_crop_limits = [-2000, 2000]; % nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGING CROPPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% Crop localizations
cropped_path = strcat(data_path,'_cropped');
figure_path = fullfile(cropped_path,'figures');

% make new path
if ~exist(cropped_path,'dir')
    mkdir(cropped_path)
end

if ~exist(figure_path,'dir')
    mkdir(figure_path)
end

files = dir(fullfile(data_path,'*.csv'));
n_files = numel(files);
file_names = {files.name}';

for i = 1:n_files

    current_file = fullfile(files(i).folder,files(i).name);
    data = readtable(current_file,'VariableNamingRule','preserve');

    if tracked
        keep_index = (z_crop_limits(1) <= data.z)...
            & (data.z <= z_crop_limits(2));
        excluded_tracks = unique(data.("#track")(~keep_index));
        data = data(~ismember(data.("#track"),excluded_tracks),:);
    elseif ~tracked
        keep_index = (z_crop_limits(1) <= data.("z (nm)"))...
            & (data.("z (nm)") <= z_crop_limits(2));
        data = data(keep_index,:);
    else
        error('User input "tracked" must be logical "ture" or "false"')
    end

    % save the data file
    writetable(data,fullfile(cropped_path,files(i).name))

    % plot data
    plotdata(data,figure_path,extractBefore(files(i).name,'.csv'),tracked)
    
end


% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotdata(data,figure_save_path,file_name,tracked) %===============
                                                                           
    if tracked
        x = data.x/1000;
        y = data.y/1000;
        z = data.z/1000;
    elseif ~tracked
        x = data.("x (nm)")/1000;
        y = data.("y (nm)")/1000;
        z = data.("z (nm)")/1000;
    end

    fig = figure(1);
    scatter3(x,y,z,3,'filled','MarkerFaceAlpha',0.2)
        
    axis equal
    axis tight

    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('z (\mum)')
      
    set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
    
    zlim([-3 3])
    set(gca,'ZTick', -3:1:3,'ZTickLabel',{'-3','','','0','','','3'})
    
    views = [ 0    90;
              0     0;
             -37.5 30];

    view_names = {'xy', 'xz','xyz'};
    
    for i = 1:numel(view_names)
        view(views(i,:))                   
        saveFigure(figure_save_path,[file_name,'_view=',view_names{i}],fig)
    end

end %======================================================================


function saveFigure(save_path,save_name,fig) %=============================
    
    full_save_name = fullfile(save_path,save_name);
    
    file_name_parts = split(full_save_name,'_view=');
    
    if strcmp(file_name_parts{2},'xy')
        savefig(fig,[file_name_parts{1},'.fig'])
    end

    saveas(fig,[full_save_name,'.png'],'png');   
    print(fig,'-vector','-dsvg',[full_save_name,'.svg'])

end %======================================================================