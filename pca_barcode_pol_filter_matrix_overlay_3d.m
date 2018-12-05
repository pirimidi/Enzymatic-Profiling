%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: April 11, 2018
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a filtered data matrix, containing all 
% unique pore IDs, which are categorized as a high-probability barcode hits
% (consensus alignment accuracy (CAA) > 80%, 50 bp < alignment length < 
% 500 bp), along with the 20 kinetic parameters (KP) to be used in the PCA 
% analysis. Next, it iterates through all unique barcodes in the measurement
% set and caluclates the PCA based on the 20 KP associated with a list of 
% pores. Finally, it generates a composite 3D plot of N PCA clusters are 
% displayed in N distinct colors, where N is the number of unique barcodes
% (which is associated with a unique POL library variant). NOTE: Use MATLAB
% R2017a to run this code, older version might trigger errors.
%
% INPUT USAGE:
%
% type = 0 - PCA is calculated for each POL-BAR
% type = 1 - the mean PCA is calculated for the POL-BAR list
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function pca_barcode_pol_filter_matrix_overlay_3d(type)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> PCA barcode-pol filter matrix start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to MAT-file directory.
mat_dir = strcat(work_dir, '\data\mat_files');
cd(mat_dir);

disp(['--> IN DIRECTORY: ' mat_dir]);

% Create empty data matrix.
DATA = [];

% Read filtered data matrix in.           
DATA = dlmread('pol6-all_data_matrix.txt');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> NORMALIZED SCATTER PLOTTING SECTION');

% Figure counter and legend array.
counter = 1; leg = {};

% Generate 'experimental set' di rectory to hold the scatter plots.
cdir = 'pca_barcode_pol_filter_matrix_overlay_3d';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end
         
% Find all barcodes (with repeats) in the data matrix.
ie = DATA(:,1);

% Filter out unique barcodes in this set.
[unique_cell_id, ia, ic] = unique(ie);

% Initialize legend counter.
j = 1;
        
% Iterate through all unique barcodes.
for i = 1:length(unique_cell_id)
    
    % Initialize data matrices.
    BAR = []; P = [];
 
    is = DATA(:,1) == unique_cell_id(i);
    BAR = DATA(is, 4:end);
    [r, c] = size(BAR);
    
    % We need at least 3 pores per barcode to calculate PCA.
    if r >=3
        
        % Increment legend.
        leg{j} = num2str(unique_cell_id(i));

        % Calculate principle components of data matrix.
        opt = statset('pca'); opt.MaxIter = 1e+6;
        [coeff, score, latent, tsquared, explained, mu] = pca(BAR, 'Algorithm', 'svd', 'Options', opt)
        P = BAR*coeff;
        
        % Check if there is only 2 principal components.
        [r, c] = size(P);
        
        % If so, expand to 3rd dimension with 0s.
        if c == 2
            P(:,3) = [0;0;0];
        end
        
        % PCA is calculated for each POL-BAR.
        if type == 0
            PC1 = P(:,1); PC2 = P(:,2); PC3 = P(:,3);
            
        % The mean PCA is calculated for the POL-BAR list.    
        else
            PC1 = mean(P(:,1)); PC2 = mean(P(:,2)); PC3 = mean(P(:,3));
        end

        % Create 3D plot of PC1 vs. PC2 vs. PC3.
        scatter3(PC1, PC2, PC3, 'o', 'MarkerEdgeColor', 'k', ...
                                     'MarkerFaceColor', rand(1,3), ...
                                     'LineWidth', 1);
        hold on;

        BAR
        size(BAR)
        explained 
        
        % Increment legend member counter.
        j = j + 1;    
    end
end
  
% Decorate plot.
grid;
title('Filtered Principal Component Analysis');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
axis([-4 10 -4 10 -4 10]);

% Insert legend.
leg
legend(leg);

% Create rotating 3D video.
OptionZ.FrameRate = 15;
OptionZ.Duration = 5.5; 
OptionZ.Periodic = true; 

% PCA is calculated for each POL-BAR.
if type == 0
    CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], 'pca_barcode_3d_rotation_filter', OptionZ);

    % The mean PCA is calculated for the POL-BAR list.    
else
    CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], 'mean_pca_barcode_3d_rotation_filter', OptionZ);
end

% Save plot.
if type == 0
    savefig('pca_barcode_pol_filter_matrix_overlay_3d.fig');
    print('-dbmp', 'pca_barcode_pol_filter_matrix_overlay_3d.bmp');
else
    savefig('mean_pca_barcode_pol_filter_matrix_overlay_3d.fig');
    print('-dbmp', 'mean_pca_barcode_pol_filter_matrix_overlay_3d.bmp');
end

% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);
movefile('*.mp4', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> PCA barcode-pol filter matrix overlay 3D end');
fprintf('\n');

% End timer.
toc