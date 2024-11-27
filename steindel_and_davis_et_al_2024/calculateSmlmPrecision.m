classdef calculateSmlmPrecision
    
    
    properties%============================================================
        
        inputPath
        
        dataPath
        dataFileNames
        dataType
        numberOfFiles
        
        precision = struct('distance2center',[],'numberOfDataPoints',[],...
                                 'numberOfLocalizationsInKernel',[])
        numberOfFrames = struct('mean',[], 'median', [],'numberOfLocalizaitons',[])      
    end%===================================================================
    
    
    methods%===============================================================
        
        
        function obj = calculateSmlmPrecision(varargin)%%%%%%%%%%%%%%%%%%%%
            
            % get input path if provided
            switch numel(varargin)
                case 0
                    obj.inputPath = '*.csv';
                case 1
                    obj.inputPath = fullfile(varargin{1},'*.csv');
                    obj.dataType = 'fixed';
                case 2
                    obj.inputPath = fullfile(varargin{1},'*.csv');
                    obj.dataType = varargin{2};
                otherwise
                    error('too many inputs! You can either input (1) nothing or (2) a data path to open up to.')
            end
            
            % Load data and set up directories
            obj = obj.getDataFileInformation;
            
            % Step through data sets and calculate precision
            disp('Stepping through data sets...')
            
            for i = 1:obj.numberOfFiles
            
                disp(['Calculating precisions for dataset: ',...
                            num2str(i),' of ', num2str(obj.numberOfFiles)])
                        
                disp('Writing precisino output for all files...')
                
                % set up saving directory for individual data file
                currentFileName = obj.dataFileNames{i};
                savePath = fullfile(obj.dataPath,...
                    strcat(extractBefore(currentFileName,'.csv'),'_precision'));
                
                if ~exist(savePath,'dir')                    
                    mkdir(savePath)                    
                end
                
                % calculate precision
                [obj_precision] = obj.calculatePrecision(currentFileName, savePath);
                
                % update parameter vectors for all files selected
                obj.precision.distance2center(i,1:size(obj_precision.distance2center.precision,2)) = obj_precision.distance2center.precision;
                obj.precision.numberOfDataPoints(i) = obj_precision.distance2center.numberOfDataPoints;
                obj.precision.numberOfLocalizationsInKernel(i) = obj_precision.distance2center.numberOfLocalizationsInKernel;
                
                obj.numberOfFrames.mean(i) = obj_precision.numberOfFrames.mean;
                obj.numberOfFrames.median(i) = obj_precision.numberOfFrames.median;
                obj.numberOfFrames.numberOfLocalizaitons(i) = obj_precision.numberOfFrames.numberOfLocalizaitons;
            end
            
            disp('...precision calculating completed')

            % write out precisions to file
            disp('Writing precisino output for all files...')
            fullSaveFileName = fullfile(obj.dataPath,'precision_summary.txt');
            
            fileID = fopen(fullSaveFileName,'w');
            
            fprintf(fileID,'\n');
            fprintf(fileID,'%15s %15s %15s %15s %15s %15s %15s %15s %50s\n',...
                            'x-precision(nm)', 'y-precision(nm)',... 
                            'z-precision(nm)','nDataPoints',... 
                            'nLocsInKernel', 'meanNumFrames',...
                            'medianNumFrames', 'nLocs',...
                            'file_name');
            fprintf(fileID,'\n');
            
            for i = 1:obj.numberOfFiles
            
                fprintf(fileID,'%15.5f %15.5f %15.5f %15.0f %15.0f %15.2f %15.2f %15.0f %50s\n',...
                    obj.precision.distance2center(i,:),...
                    obj.precision.numberOfDataPoints(i),...
                    obj.precision.numberOfLocalizationsInKernel(i),...
                    obj.numberOfFrames.mean(i),...                    
                    obj.numberOfFrames.median(i),...
                    obj.numberOfFrames.numberOfLocalizaitons(i),...
                    obj.dataFileNames{i});                            
                
            end
            
            fclose(fileID);
            
            disp('...Output written.')
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        function obj = getDataFileInformation(obj)%%%%%%%%%%%%%%%%%%%%%%%%%
                                   
            % get location of localizations using dialog box---------------            
            disp('User input initiated...')
            
            disp('Please select a localizations file...')
            
            % getting selected files
            [localizations_file, localizations_path] = ...
                uigetfile(obj.inputPath,'Select localizations file',...
                                                       'MultiSelect','on');            
            
            % make sure file names are of tyoe cell
            if ~iscell(localizations_file)
                localizations_file = {localizations_file};
            end
            
            % updat object properties
            obj.dataFileNames = localizations_file;                                                    
            obj.dataPath = localizations_path;                                    
            obj.numberOfFiles = numel(obj.dataFileNames);            
            obj.precision.distance2center = zeros(obj.numberOfFiles,3);
            obj.precision.numberOfDataPoints = zeros(obj.numberOfFiles,1);                       
            obj.precision.numberOfLocalizationsInKernel = zeros(obj.numberOfFiles,1);
            obj.numberOfFrames.mean = zeros(obj.numberOfFiles,1);
            obj.numberOfFrames.median = zeros(obj.numberOfFiles,1);
            obj.numberOfFrames.numberOfLocalizaitons = zeros(obj.numberOfFiles,1);
            
            disp('...user input completed')
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
        function obj_precision = calculatePrecision(obj,currentFileName,savePath)
                
            % Get data for current file
            [localization_id, coordinates] = loadData(obj,currentFileName);
            
            switch obj.dataType
                case 'fixed'
                % perform precision algorithm
                obj_precision = precisionCalculator(localization_id, coordinates,savePath);                                                                                             

                % plot kernels            
                obj_precision.plotDistance2center();
                obj_precision.plotNumberOfFrames();
                case 'live'
                otherwise
                    
            end
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
          
        function [localization_id, coordinates] = loadData(obj,currentFileName)
                                  
            file = importdata(fullfile(obj.dataPath,currentFileName));
            
            localization_id = file.data(:,1);
            coordinates = file.data(:,3:end);
                                  
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    end%===================================================================
    
    
end%=======================================================================
%==========================================================================
