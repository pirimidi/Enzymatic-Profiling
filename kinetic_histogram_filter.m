%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: April 7, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters parsed from the 'cell_annotations' 
% EXCEL file. Then, for each 'pol-bar', it iterates through the experimental 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates normalized 
% histograms for each filtered 'pol-bar' case displaying them on the same 
% plot for comparison. NOTE: Use MATLAB R2017a to run this code, older 
% version might trigger errors. Vary 'pol', fix 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_histogram_filter(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic histogram filter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

% Color palette.
col = hsv(3); col_id = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
data_dir = strcat(work_dir, '/data/mat_files/merge');
cd(data_dir);

disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_merge*.mat');

% Figure counter.
counter = 1;

disp('--> HISTOGRAM PLOTTING SECTION');
    
% Iterate through all three barcodes.
for bar = 1:length(barcodes)
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    %                            HISTOGRAM PLOTTING                           %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    % Generate 'experimental set' directory to hold the whisker plots.
    cdir = ['kinetic_histogram_filter_', barcodes{bar}];
    
    % Define direcory to hold figures.
    if ~exist(cdir, 'dir')
      mkdir(cdir);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TTT median histogram.
    coi_name = 'ttt_median';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                    
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
            
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered TTT Histogram');
        xlabel('Time-to-Thread (ms)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 4 0 0.5]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_ttt_median_filtered'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean dwell time histogram for the A-tagged nucleotide.
    coi_name = 'level_call_A_mean_dwell_time';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered Mean Dwell Time Histogram');
        xlabel('Dwell Time (s)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 6 0 1]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_A_dwell_time_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean dwell time histogram for the C-tagged nucleotide.
    coi_name = 'level_call_C_mean_dwell_time';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered Mean Dwell Time Histogram');
        xlabel('Dwell Time (s)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 6 0 1]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_C_dwell_time_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean dwell time histogram for the G-tagged nucleotide.
    coi_name = 'level_call_G_mean_dwell_time';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered Mean Dwell Time Histogram');
        xlabel('Dwell Time (s)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 6 0 1]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_G_dwell_time_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean dwell time histogram for the T-tagged nucleotide.
    coi_name = 'level_call_T_mean_dwell_time';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered Mean Dwell Time Histogram');
        xlabel('Dwell Time (s)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 6 0 1]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_T_dwell_time_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Median TTT histogram for the A-tagged nucleotide.
    coi_name = 'level_call_A_ttt_median';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered A-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 4 0 0.45]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_A_ttt_median_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Median TTT histogram for the C-tagged nucleotide.
    coi_name = 'level_call_C_ttt_median';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered C-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 4 0 0.45]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_C_ttt_median_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Median TTT histogram for the G-tagged nucleotide.
    coi_name = 'level_call_G_ttt_median';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered G-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 4 0 0.45]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_G_ttt_median_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Median TTT histogram for the T-tagged nucleotide.
    coi_name = 'level_call_T_ttt_median';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1

                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                ci = extractfield(pb_f3, coi_name);

                % Generate container arrays for all experiments.
                COI = [COI; ci'];
                
            end
        end

        % Determine number of total events in this data set.
        ll = ['N = ' , num2str(length(COI))];

        % Decorate accumulative histogram.
        histogram(COI, 10, 'FaceColor', col(col_id,:), 'Normalization', 'probability');
        grid;
        title('Filtered T-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
        ylabel('Probability');
        legend(ll, 'Location', 'eastoutside');
        axis([0 4 0 0.45]);

        % Save plot.
        ss = strsplit(list(e).name, '.');
        pp = strsplit(ss{1}, '_');
        fn = [num2str(counter), '_', pp{1}, '_', barcodes{bar}, '_T_ttt_median_filter'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;
            
        % Increment color ID.
        col_id = col_id + 1;

        % Reset color ID.
        if col_id == 4
            col_id = 1;
        end    

        disp(['--> PROCESSED DATA SET: ' fn]); 

    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Move all figures to 'plots' directory.
    movefile('*.fig', cdir);
    movefile('*.bmp', cdir);

    % Close all open figures.
    close all;
    
end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Kinetic alignment filter end');
fprintf('\n');

% End timer.
toc