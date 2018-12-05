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

function pca_mapback_pol_filter_matrix_overlay_3d(num_sets)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> PCA map back filter matrix start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

col = {'b', 'k', 'r'};
col_id = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to MAT-file directory.
mat_dir = strcat(work_dir, '\data\final_set\mat_files');
cd(mat_dir);

disp(['--> IN DIRECTORY: ' mat_dir]);

% Read in a TXT file for predefined 96 naive barcodes.
bar_seqs = fastaread('unique_barcodes_50%.txt');

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_*merge*.mat');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> NORMALIZED SCATTER PLOTTING SECTION');

% Figure counter.
counter = 1;

% Create empty data matrix.
DATA = [];

% Generate 'experimental set' directory to hold the scatter plots.
cdir = 'pca_mapback_pol_filter_matrix_overlay_3d';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Read filtered data matrix in.           
DATA = dlmread('data_matrix.txt');
            
for s = 1:num_sets 

    % Clear container.
    %BAR = [];

    % Analysis type 1 - barcodes 1-32 considered in PCA analysis.
    if s == 1

        % All barcodes in set 1.
        i1 = DATA(:,1)<33;
        BAR = DATA(i1, 4:end);

    % Analysis type 2 - barcodes 33-64 considered in PCA analysis.
    elseif s == 2

        % All barcodes in set 2.
        ii = DATA(:,1)>33;
        DATA2 = DATA(ii,:);
        i2 = DATA2(:,1)<65;
        BAR = DATA2(i2, 4:end);

    % Analysis type 3 - barcodes 65-96 considered in PCA analysis.
    elseif s == 3

        % All barcodes in set 3.
        i3 = DATA(:,1)>64;
        BAR = DATA(i3, 4:end);
    
    end

    % Calculate principle components of data matrix.
    opt = statset('pca'); opt.MaxIter = 1e+6;
    [coeff, score, latent, tsquared, explained, mu] = pca(BAR, 'Algorithm', 'svd', 'Options', opt);
    P = BAR*coeff;
    PC1 = P(:,1); PC2 = P(:,2); PC3 = P(:,3);

    % Create 3D plot of PC1 vs. PC2 vs. PC3.
    scatter3(PC1, PC2, PC3, 'o', 'MarkerEdgeColor', 'k', ...
                                 'MarkerFaceColor', col{col_id}, ...
                                 'LineWidth', 1);
    hold on;

    % Increment color ID.
    col_id = col_id + 1;

    size(BAR)
    explained

    disp('--> PROCESSED DATA SET: data_matrix.txt'); 
        
 end    
  
% Decorate plot.
grid;
title('Filtered Principal Component Analysis');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
axis([-4 10 -4 10 -4 10]);

% Insert legend.
legend('CBT01-32', 'CBT33-64', 'CBT65-96');

% Create rotating 3D video.
OptionZ.FrameRate = 15;
OptionZ.Duration = 5.5; 
OptionZ.Periodic = true; 
CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], 'pca_mapback_3d_rotation_filter', OptionZ);

% Save plot.
savefig('pca_mapback_pol_filter_matrix_overlay_3d.fig');
print('-dbmp', 'pca_mapback_pol_filter_matrix_overlay_3d.bmp');

% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);
movefile('*.mp4', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> PCA mapback pol filter matrix overlay 3D end');
fprintf('\n');

% End timer.
toc