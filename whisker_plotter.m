%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: February 15, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters parsed from the 'cell_annotations' 
% EXCEL file. Then, for each 'pol-bar', it iterates through the experimental 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates whisker plots 
% for each experiment (components of total) and a final, cumulative whisker
% plot all experiments (all data). NOTE: Use MATLAB R2017a to run this 
% code, older version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function whisker_plotter

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Whisker plotter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

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

disp('--> WHISKER PLOTTING SECTION');
    
% Iterate through all 'experiment' folders.
for e = 1:length(list) 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    %                              WHISKER PLOTTING                           %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % Load nanopore data structure per experimental set (pol-bar).
    load(list(e).name);
    
    fprintf('\n');
    disp(['--> MAT-FILE LOADED: ' list(e).name]);
    
    % Generate 'experimental set' directory to hold the whisker plots.
    ss = strsplit(list(e).name, '.');
    cdir = ['whisker_plot_', ss{1}];
    
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
            
            % Define array container for 'coi'.
            COI = [];

            % Generate figure for each pol-bar combination.
            figure(counter);

            % Iterate through all 'experiment' folders for (pol-)bar.
            for exp = 1:length(cells_list(bar).experiment)

                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                values = extractfield(pb, coi_name);

                % Total number of single pores and its range.
                seq_pores = length(pb); x = 1:seq_pores;

                % Generate container arrays for all experiments.
                COI = [COI; values'];
                
                % Determine number of events in this data set.
                SE(exp) = length(values);
                
            end

            % Update counter.
            counter = counter + 1;

            % Build grouping variable array.
            GE = []; 

            for j = 1:length(cells_list(bar).experiment)

                % Initialize containers.
                ge = [];

                for k = 1:SE(j)
                    ge{k} = j;
                end

                % Update concatinated grouping variable array.
                GE = [GE; ge'];

            end

            % Create dwell time (s) boxplot of all filtered events (exp).
            boxplot(COI, GE, 'PlotStyle', 'compact'); 
            grid;
            title([ordinate_name, ' whisker plot']);
            xlabel('cell id');
            ylabel(ordinate_name);

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
end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> whisker plotter end');
fprintf('\n');

% End timer.
toc