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
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates whisker plots 
% featuring kinetic differences for each filtered 'pol-bar' case displaying
% them on the same plot for comparison. NOTE: Use MATLAB R2017a to run this
% code, older version might trigger errors. Fix 'pol', vary 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_whisker_pol_filter(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic whisker pol filter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

% (6) No upper/lower bound 'coi'.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)'}; 
 
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

disp('--> WHISKER PLOTTING SECTION');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              WHISKER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate 'experimental set' directory to hold the whisker plots.
cdir = 'kinetic_whisker_pol_filter';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Iterate through all 'cois'.
for coi_id = 1:length(coi)

    % Ordinate name.
    coi_split = strsplit(coi{coi_id});
    coi_name = coi_split{1}
    ordinate_name = regexprep(coi_name, '_', '-');
    
    % Define array container for fixed 'bar' and grouping variable array.
    POL = []; GE = []; 
        
    % Iterate through all 'pol' folders.
    for e = 1:length(list) 
        
        % Define array container for 'coi' and grouping variable array.
        COI = []; ge = []; 
           
        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        ss = strsplit(list(e).name, '.');
            
        % Iterate through all three barcodes.
        for bar = 1:length(barcodes)
            
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
   
            % Generate container array for 'pol-bar' combination (fix bar).
            POL = [POL; COI]; 
            
            % Update concatinated grouping variable array.
            for k = 1:length(COI)
                ge{k} = bar;
            end

            GE = [GE; ge'];
                 
        end
    end
    
    % Generate figure for each pol-bar combination.
    figure(counter);

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(POL, GE, 'PlotStyle', 'compact', 'Symbol', 'b', 'Whisker', 0.2413);
    grid;
    title([ordinate_name, ' pol-kinetics']);
    xlabel('bar id');
    ylabel(ordinate_name);

    % Set y-axis limits.
    if strcmp(coi_name, 'level_call_A_k_cat_rate') || strcmp(coi_name, 'level_call_C_k_cat_rate') || ...
       strcmp(coi_name, 'level_call_T_k_cat_rate') || strcmp(coi_name, 'level_call_G_k_cat_rate')

        ylim([0 4.5]);

    end

    if strcmp(coi_name, 'level_call_A_k_off_rate') || strcmp(coi_name, 'level_call_C_k_off_rate') || ...
       strcmp(coi_name, 'level_call_T_k_off_rate') || strcmp(coi_name, 'level_call_G_k_off_rate')

        ylim([-1.5 1]);

    end

    if strcmp(coi_name, 'level_call_A_mean_dwell_time') || strcmp(coi_name, 'level_call_C_mean_dwell_time') || ...
       strcmp(coi_name, 'level_call_T_mean_dwell_time') || strcmp(coi_name, 'level_call_G_mean_dwell_time')

        ylim([0 2.5]);

    end

    % Save plot.
    fn = [num2str(coi_id), '_', coi_name];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);

    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);   
     
end
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;
    
% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Kinetic whisker pol filter end');
fprintf('\n');

% End timer.
toc