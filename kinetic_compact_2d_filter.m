%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: July 27, 2017
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
% code, older version might trigger errors. Vary 'pol', fix 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_compact_2d_filter(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic compact 2D filter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

% Dwell time and k_cat parameters COIs for each of the four bases.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_mean_dwell_time (seconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_mean_dwell_time (seconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_mean_dwell_time (seconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_mean_dwell_time (seconds)'}; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
data_dir = strcat(work_dir, '/data/mat_files/merge');
cd(data_dir);
DATA
disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_merge*.mat');

disp('--> COMPACT 2D PLOTTING SECTION');
    
% Iterate through all three barcodes.
for bar = 1:length(barcodes)
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    %                           COMPACT 2D PLOTTING                           %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    % Generate 'experimental set' directory to hold the whisker plots.
    cdir = ['kinetic_compact_2d_filter_', barcodes{bar}];
    
    % Define direcory to hold figures.
    if ~exist(cdir, 'dir')
      mkdir(cdir);
    end

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 
        
        % Define array container for fixed 'bar'.
        BAR = [];

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
            
        % Iterate through all 'cois'.
        for coi_id = 1:length(coi)
        
            % Ordinate name.
            coi_split = strsplit(coi{coi_id});
            coi_name = coi_split{1}
            ordinate_name = regexprep(coi_name, '_', '-');
    
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

            % Generate container array for 'pol-bar' combination (fix bar).
            BAR = horzcat(BAR, COI);
            
        end     
        
        % Write COI data to a matrix and save it.           
        dlmwrite(['pol', num2str(e), '_', barcodes{bar}, '_matrix.txt'], BAR, 'delimiter', '\t', 'precision', '%.4f');
        size(BAR)
        BAR 
        
        disp(['--> PROCESSED DATA SET: ' list(e).name]);     
    end
    
    % Move all text files to directory.
    movefile('*_matrix.txt', cdir);
   
end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Kinetic compact 2D filter end');
fprintf('\n');

% End timer.
toc