%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: April 5, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all alignment parameters generated by 
% 'excel_parser.m'. Then, for each 'pol-bar', it iterates through 
% the experimental structures and generates a cumulative matrix of nanopore
% measurements (row) vs. characteristic of interests (column). Subsequently,
% it calculates the PCA coefficient matrix, using the default settings,
% which is then used to compute PC1 and PC2 set of values. Finally, 2D 
% scatter plots are generated for each 'pol-bar' case to feature unique 
% grouping. NOTE: Use MATLAB R2017a to run this code, older version might 
% trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function pca_kinetic_pol_filter(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> PCA kinetic pol filter start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

col = {'r', 'b', 'k', 'y', 'm', 'c', 'g', [255,128,0]./255, [127,0,255]./255};
col_id = 1;

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
list = dir('pol6*_*merge*.mat');

% Figure counter.
counter = 1;

disp('--> NORMALIZED SCATTER PLOTTING SECTION');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate 'experimental set' directory to hold the scatter plots.
cdir = 'pca_kinetic_pol_filter';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Iterate through all 'pol' folders.
 for e = 1:length(list) 
     
    % Define array container for fixed 'bar' array.
    BAR = [];

    % Iterate through all three barcodes.
    for bar = 1:length(barcodes)  

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        ss = strsplit(list(e).name, '.');

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = merge_list(bar).experiment(exp).s;
            
            % Define array container for 'coi' array.
            COI = [];
            
            % Iterate through all 'coi'.
            for coi_id = 1:length(coi)

                % Ordinate name.
                coi_split = strsplit(coi{coi_id});
                coi_name = coi_split{1};
                
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
                COI = horzcat(COI, ci');
                
            end 
                        
            % Generate container array for 'pol-bar' combination (fix bar).
            BAR = vertcat(BAR, COI);
        
        end        
    end 
    
    % Calculate principle components of data matrix.
    opt = statset('pca'); opt.MaxIter = 1e+6;
    [coeff, score, latent, tsquared, explained, mu] = pca(BAR, 'Algorithm', 'svd', 'Options', opt);
    P = BAR*coeff;
    PC1 = P(:,1); PC2 = P(:,2);

    % Generate figure for each pol-bar combination.
    figure(counter);

    % Create 2D plot of PC1 vs. PC2.
    scatter(PC1, PC2, 'o', 'MarkerEdgeColor', 'k', ...
                           'MarkerFaceColor', col{col_id}, ...
                           'LineWidth', 1);
    grid;
    title('Filtered Principal Component Analysis');
    xlabel('PC1');
    ylabel('PC2');
    axis([-4 8 -4 10]);

    % Generate name/legend for figure.
    pp = strsplit(ss{1}, '_');
    fn = [pp{1}, '_pca_filter'];
    legend(pp{1});

    % Save plot.
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);

    % Update counter.
    counter = counter + 1;

    % Increment color ID.
    col_id = col_id + 1;

    % Reset color ID.
    if col_id == 21
        col_id = 1;
    end

    size(BAR)
    explained

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
disp('--> PCA kinetic pol filter end');
fprintf('\n');

% End timer.
toc