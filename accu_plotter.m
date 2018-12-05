%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: February 10, 2017 > my birthday! :)
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters parsed from the 'cell_annotations' 
% EXCEL file. Then, for each 'pol-bar', it iterates through the experimental 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates accumulative 
% 2D plots overlaying each data set on the same plot. NOTE: Use MATLAB 
% R2017a to run this code, older version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function accu_plotter

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Accumulative plotter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};
%barcodes = {'fv2'};

% Color palette.
col = hsv(20); col_id = 1;

% Create an pre-populated cell array for keyset in hash table.
coi = {'row (#)', 'col (#)', ...
       'sequencing_end_time (seconds)', 'sequencing_lifetime (seconds)', 'lifetime_after_tag_flow (seconds)', 'len_starting_single_pore_run (#)', ...
       'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'align_copies (#)', ...
       'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
       'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)', ...
       'level_call_A_counts (#)', 'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_counts (#)', 'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_counts (#)', 'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_counts (#)', 'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
data_dir = strcat(work_dir, '\data\mat_files');
cd(data_dir);

disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_stats.mat');

% Figure counter.
counter = 1;

disp('--> ACCUMULATIVE PARAMETER PLOTTING SECTION');
    
% Iterate through all 'experiment' folders.
for e = 1:length(list) 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    %                    ACCUMULATIVE PARAMETER PLOTTING                      %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % Load nanopore data structure per experimental set (pol-bar).
    load(list(e).name);
    
    fprintf('\n');
    disp(['--> MAT-FILE LOADED: ' list(e).name]);
    
    % Generate 'experimental set' directory to hold the accumulative plots.
    ss = strsplit(list(e).name, '.');
    cdir = ['plots_', ss{1}];
    
    % Define direcory to hold figures.
    if ~exist(cdir, 'dir')
      mkdir(cdir);
    end

    % Iterate through all 'cois'.
    for coi_id = 1:length(coi)
        
        % Ordinate name.
        coi_split = strsplit(coi{coi_id});
        coi_name = coi_split{1};
        ordinate_name = regexprep(coi_name, '_', '-');
                
        % Iterate through all three barcodes.
        for bar = 1:length(cells_list)

            % Generate figure for each pol-bar combination.
            figure(counter); ll = {};

            % Iterate through all 'experiment' folders for (pol-)bar.
            for exp = 1:length(cells_list(bar).experiment)

                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Sort cells as desired for 'coi'.
                cells_sorted = nestedSortStruct2(pb, coi_name, -1);
                values = extractfield(cells_sorted, coi_name);
                keys = extractfield(cells_sorted, 'cell_id');

                % Total number of single pores and its range.
                seq_pores = length(pb); x = 1:seq_pores;

                % Barcode iteration vs. cell id plot.
                plot(x, values, 'rx', 'MarkerSize', 2, 'color', col(col_id,:));
                hold on;

                % Increment legend label.
                ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
                
                % Increment color ID.
                col_id = col_id + 1;

                % Reset color ID.
                if col_id == 21
                    col_id = 1;
                end
            end

            hold off;

            % Update counter.
            counter = counter + 1;

            % Decorate overlaid experimental plot.	
            grid;
            title([ordinate_name, ' vs cell id']);
            xlabel('cell id');
            ylabel(ordinate_name);
            legend(ll, 'Location', 'eastoutside');

            % Save plot.
            pp = strsplit(ss{1}, '_');
            fn = [num2str(coi_id), '_', pp{1}, '_', barcodes{bar}, '_', coi_name];
            savefig([fn, '.fig']);
            print('-dbmp', [fn, '.bmp']);

            disp(['--> PROCESSED DATA SET: ' fn]);     
        end

        % Move all figures to 'plots' directory.
        movefile('*.fig', cdir);
        movefile('*.bmp', cdir);

        % Close all open figures.
        close all;
    end

    % Read length histogram.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'align_read_length');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % Read length histogram.
            histogram(values, length(values), 'DisplayStyle', 'stairs', 'EdgeColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('Read Length Histogram');
        xlabel('Read Length (bp)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+1), '_', pp{1}, '_', barcodes{bar}, '_read_length_hist'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end
    
    % Single pore lifetime histogram.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'sequencing_lifetime');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % Sequencing lifetime.
            histogram(values, 20, 'FaceColor', col(col_id,:));
            xlim([0 7200]);
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('Single Pore Lifetime Histogram');
        xlabel('Single Pore Lifetime (s)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
        %set(gca, 'YScale', 'log');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+2), '_', pp{1}, '_', barcodes{bar}, '_sequencing_lifetime'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end

    % TTT histogram for the 4 tagged nucleotides.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'ttt_median');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % All tags.
            histogram(values, 10, 'FaceColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('Accumulative TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+3), '_', pp{1}, '_', barcodes{bar}, '_ttt_median'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end

    % TTT histogram for the C-tagged nucleotide.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'level_call_C_ttt_median');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % C-tag.
            histogram(values, 10, 'FaceColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('C-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+4), '_', pp{1}, '_', barcodes{bar}, '_level_call_C_ttt_median'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end
    
    % TTT histogram for the A-tagged nucleotide.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'level_call_A_ttt_median');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % A-tag.
            histogram(values, 10, 'FaceColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('A-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+5), '_', pp{1}, '_', barcodes{bar}, '_level_call_A_ttt_median'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end
    
    % TTT histogram for the T-tagged nucleotide.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'level_call_T_ttt_median');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % C-tag.
            histogram(values, 10, 'FaceColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('T-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+6), '_', pp{1}, '_', barcodes{bar}, '_level_call_T_ttt_median'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end
    
    % TTT histogram for the G-tagged nucleotide.
    for bar = 1:length(cells_list)

        % Generate figure for each pol-bar combination.
        figure(counter); ll = {};

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;
            values = extractfield(pb, 'level_call_G_ttt_median');

            % Total number of single pores and its range.
            seq_pores = length(pb); x = 1:seq_pores;

            % C-tag.
            histogram(values, 10, 'FaceColor', col(col_id,:));
            hold on;

            % Increment legend label.
            ll{exp} = ['data set ', num2str(exp), ': N = ', num2str(seq_pores)];
            
            % Increment color ID.
            col_id = col_id + 1;

            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end

        hold off;

        % Update counter.
        counter = counter + 1;

        % Decorate overlaid experimental plot.  
        grid;
        title('G-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend(ll, 'Location', 'eastoutside');
       
        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [num2str(length(coi)+7), '_', pp{1}, '_', barcodes{bar}, '_level_call_G_ttt_median'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        disp(['--> PROCESSED DATA SET: ' fn]);     
    end

    % Move all figures to 'plots' directory.
    movefile('*.fig', cdir);
    movefile('*.bmp', cdir);

    % Close all open figures.
    close all;
end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Accumulative plotter end');
fprintf('\n');

% End timer.
toc