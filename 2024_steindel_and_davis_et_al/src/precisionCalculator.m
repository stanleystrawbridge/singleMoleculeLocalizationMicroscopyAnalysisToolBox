classdef precisionCalculator

    
    properties%============================================================
        
        % Raw data
        localizaitonIdVector        
        coordinates
                
        uniqueIds % labels
        index     % indices of labes
        
        axLimits  % data rounded up data extrema
        
        savePath        
        
        % precision measurement
        distance2center = struct('precision',[],'numberOfDataPoints',[],...
                                 'numberOfLocalizationsInKernel',[],'kernel',[]);
        
        % parameters of the particle tracks through time
        numberOfFrames = struct('mean',[], 'median', [],...
                                'numberOfLocalizaitons',[], 'kernel', [])
                
    end%===================================================================
    
    
    methods%===============================================================
        
        
        function obj_precision =  precisionCalculator(localization_id,coordinates, varargin)
            
            % updat object properties
            obj_precision.localizaitonIdVector = localization_id;
            [obj_precision.uniqueIds, ~, obj_precision.index] = ...
                                unique(obj_precision.localizaitonIdVector);
            obj_precision.coordinates = coordinates;
             
            % get Save path if provided
            switch numel(varargin)
                case 0
                    obj_precision.savePath = pwd;
                case 1
                    obj_precision.savePath =varargin{1};
                otherwise
                    error('too many inputs! You can either input (1) nothing or (2) a save path.')
            end
            
            % calculate precision by distance to center of mass of tracked
            % localization across frames
            obj_precision = obj_precision.precisionDistance2center;
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function obj_precision = precisionDistance2center(obj_precision)%%%%
                       
            % Get center of mass for each particle
            [mu, numberOfFrames] = ...
                grpstats(obj_precision.coordinates,obj_precision.index, {'mean', 'numel'});                        
            
            % save fram number parameters
            obj_precision.numberOfFrames.mean = mean(numberOfFrames(:,1));
            obj_precision.numberOfFrames.median = median(numberOfFrames(:,1));
            obj_precision.numberOfFrames.kernel = numberOfFrames(:,1);
            obj_precision.numberOfFrames.numberOfLocalizaitons = size(mu,1);
            
            % move cm of all tracks to origin
            translated_coordinates = obj_precision.coordinates - ...
                                               mu(obj_precision.index,:);
            
            % remove all singleton trajectories
            singletons = all(translated_coordinates == 0, 2);            
            translated_coordinates(singletons,:) = [];                       
            
            numberOfDataPoints = sum(~singletons);
            numberOfLocalizationsInKernel = numel(unique(obj_precision.index(~singletons)));
                        
            % Store precision kernel  
            obj_precision.distance2center.kernel = translated_coordinates;
            obj_precision.distance2center.numberOfDataPoints = numberOfDataPoints;
            obj_precision.distance2center.numberOfLocalizationsInKernel = numberOfLocalizationsInKernel;
            
            % Fite normals to idivual axes (not interested in covariance at
            % the moment), could have also just taken standard deviation...
            [muHat,sigmaHat, muCI, sigmaCI] = normfit(translated_coordinates);
            
            % Store precission: the standard deviation of the 
            obj_precision.distance2center.precision = sigmaHat;
            
            % Set up axix limits for plotting, rounding up to the closes 50
            [lower_bounds, upper_bounds] = bounds(obj_precision.distance2center.kernel);     
                                    
            obj_precision.axLimits = ceil(max(abs([lower_bounds;upper_bounds])/50))*50;                       
            obj_precision.axLimits = [ -ceil(max(abs([lower_bounds;upper_bounds])/50))*50;
                                        ceil(max(abs([lower_bounds;upper_bounds])/50))*50];
            obj_precision.axLimits = reshape(obj_precision.axLimits,1,6);
                        
            % Write precision fitting to file
            fullSaveFileName = fullfile(obj_precision.savePath,'precision_summary.txt');
            fileID = fopen(fullSaveFileName,'w');
            
            fprintf(fileID,'%20s %20s\n', 'numberOfDatapoints', 'numberOfLocalizationsInKernel');
            fprintf(fileID,'%20.0f %20.0f\n',...
                obj_precision.distance2center.numberOfDataPoints,...
                obj_precision.distance2center.numberOfLocalizationsInKernel);
            fprintf(fileID,'\n');

            fprintf(fileID,'Mean (nm)\n');
            fprintf(fileID,'%8.5f %8.5f %8.5f\n', muHat);
            fprintf(fileID,'\n');
            
            fprintf(fileID,'Mean Confidence Interval (nm)\n');
            fprintf(fileID,'%8.5f %8.5f %8.5f\n', muCI);
            fprintf(fileID,'\n');
            
            fprintf(fileID,'\n');
            fprintf(fileID,'Standard Deviation (nm)\n');
            fprintf(fileID,'%8.5f %8.5f %8.5f\n', sigmaHat);
            fprintf(fileID,'\n');
            
            fprintf(fileID,'Standard Deviation Confidence Interval (nm)\n');
            fprintf(fileID,'%8.5f %8.5f %8.5f\n', sigmaCI);
            fprintf(fileID,'\n');
            
            fprintf(fileID,'Number of Frames particle was observed\n');
            fprintf(fileID,'%20s %20s %20s\n', 'mean', 'median', 'numberOfLocalizatoins');           
            fprintf(fileID,'%20.2f %20.2f %20.0f\n',...
                                    obj_precision.numberOfFrames.mean,...
                                    obj_precision.numberOfFrames.median,...
                                    obj_precision.numberOfFrames.numberOfLocalizaitons);           
            
            fclose(fileID);
            
            % Write kernels to csv
            fullSaveFileName = fullfile(obj_precision.savePath,'kernel_distances.txt');
            writematrix(obj_precision.distance2center.kernel,fullSaveFileName)
            
            fullSaveFileName = fullfile(obj_precision.savePath,'kernel_numberOfFrames.txt');
            writematrix(obj_precision.numberOfFrames.kernel, fullSaveFileName)
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
               
        function plotDistance2center(obj_precision)%%%%%%%%%%%%%%%%%%%%%%%%                       
            
            % Plot 3D scatter plot of distances----------------------------
            fig = figure;
            clf
           
            % 3D-scatter plot of precision kernel
            scatterDistance2center(obj_precision)
            view([50 10]) % xyz-view
            title('Distance to Center')
            savefigure(obj_precision, fig, 'scatter_3d')
            
            % xy-scatter plot of precision kernel
            title('Distance to Center')
            view([0 90]) % xy-view
            savefigure(obj_precision, fig, 'scatter_xy-view')

            % xz-scatter plot of precision kernel
            title('Distance to Center')
            view([0 0]) % xz-view
            savefigure(obj_precision, fig, 'scatter_xz-view')
            
            % yz-scatter plot of precision kernel
            title('Distance to Center')
            view([90 0]) % yz-view
            savefigure(obj_precision, fig, 'scatter_yz-view')
            
            % Plot Marginal distributions of the distances-----------------
            
            % Univariate histogram of x-distances
            clf
            histogramDistance2center(obj_precision,1)   
            title('x-Distance to Center')
            savefigure(obj_precision, fig, 'histogram_x-precision')
            
            % Univariate histogram of y-distances
            clf
            title('y-Distance to Center')
            histogramDistance2center(obj_precision,2)              
            savefigure(obj_precision, fig, 'histogram_y-precision')
            
            % Univariate histogram of z-distances
            clf
            title('z-Distance to Center')
            histogramDistance2center(obj_precision,3)              
            savefigure(obj_precision, fig, 'histogram_z-precision')
            
            close all hidden
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function scatterDistance2center(obj_precision)%%%%%%%%%%%%%%%%%%%%%
                                            
            % make xyz-axes as thick line
            plot3(obj_precision.axLimits(1:2), [0 0], [0 0],'k',...
                                                            'LineWidth',2)
            hold on                        
            plot3([0 0], obj_precision.axLimits(3:4), [0 0],'k',...
                                                            'LineWidth',2)
            plot3([0 0], [0 0], obj_precision.axLimits(5:6),'k',...
                                                            'LineWidth',2)                       
            
            % scater all the data points
            scatter3(obj_precision.distance2center.kernel(:,1),...
                     obj_precision.distance2center.kernel(:,2),...
                     obj_precision.distance2center.kernel(:,3),...
                     25,'k','filled','MarkerFaceAlpha',0.1)
            
            % plot the mean
            scatter3(mean(obj_precision.distance2center.kernel(:,1)),...
                     mean(obj_precision.distance2center.kernel(:,2)),...
                     mean(obj_precision.distance2center.kernel(:,3)),...
                     75,'r','filled')
                 
            % Make it look pretty
            axis(obj_precision.axLimits) 
            axis equal
            grid on
            
            xlabel('x-distance')
            ylabel('y-distance')            
            zlabel('z-distance')  
            
            set(gca,'FontSize',12,'FontWeight','bold')
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function histogramDistance2center(obj_precision,dim)%%%%%%%%%%%%%%%%%%%%%
            
            % get univariate data
            data = obj_precision.distance2center.kernel(:,dim);
                                
            % Plot univariate data                    
            h = histogram(data,'FaceColor','k', 'FaceAlpha',0.5,...
                               'DisplayStyle','bar','EdgeColor','none',...
                               'Normalization','pdf');                                   
            hold on  
             
            % Get the axis limits
            idx = 2*(dim - 1)+1;           
            x_limits = obj_precision.axLimits(idx:idx+1);                 
                         
            % Generate univariate normal distribution from normal
            % parameters
            
            % Calculate parameters
            mu = mean(data); 
            sigma = std(data);
            
            % generate normal pdf
            x_normal = linspace(x_limits(1), x_limits(2),100);
            y_normal = normpdf(x_normal,mu,sigma);
                                 
            % Plot normal distribution
            plot(x_normal,y_normal,'k', 'LineWidth',3)            
            
            % Plot mean line
            yMax = max([h.Values, y_normal]);
            plot( mu*[1 1], [0 yMax], 'r--', 'LineWidth', 2)           
            
            % Make it look pretty
            xlim(x_limits);
            ylim([0, yMax])           
            grid on
            
            xlabel('distance')
            ylabel('count')
            
            set(gca,'FontSize',12,'FontWeight','bold')
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function plotNumberOfFrames(obj_precision)%%%%%%%%%%%%%%%%%%%%%%%%
            
            fig = figure;
            clf
           
            % Plot univariate data
            h = histogram(obj_precision.numberOfFrames.kernel,...
                        'FaceColor','k', 'FaceAlpha',0.5,...
                        'DisplayStyle','bar','EdgeColor','none',...
                        'Normalization','pdf');
            hold on
            
            % Get the axis limits
            x_limits = [0, max(obj_precision.numberOfFrames.kernel)];
            y_max = max(h.Values);
            y_limits = [0 y_max];
            
            % Plot mean and median lines
            plot(obj_precision.numberOfFrames.mean*[1 1],[0 y_max],...
                                                        'r--','LineWidth',2)             
            plot(obj_precision.numberOfFrames.median*[1 1],[0 y_max],...
                                                    'k--','LineWidth',2) 
             
            legend('kernel','mean','median')
            % Make it look pretty
            axis([x_limits y_limits]);              
            grid on
            
            title('Tracked Localizations')
            xlabel('Number of Frames Observed')
            ylabel('count')
            
            set(gca,'FontSize',12,'FontWeight','bold')
            
            savefigure(obj_precision, fig, 'histogram_numberOfFrames')
            
            close all hidden
            
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function savefigure(obj_precision, fig, saveName)
            
            full_save_name = fullfile(obj_precision.savePath, saveName);
                        
            saveas(fig, full_save_name,'fig')
            saveas(fig, full_save_name,'tif')
            saveas(fig, full_save_name,'svg')        
   
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
        
    end%===================================================================
    
    
end%=======================================================================
%==========================================================================


