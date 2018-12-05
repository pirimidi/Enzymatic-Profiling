%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: Jult 24, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a MAT-file, containing all alignment 
% parameters generated by 'consensus_alignment_96.m'. Then, it iterates 
% through the experimental structures and collects all pore IDs, which are 
% categorized as a high-probability barcode hits (consensus alignment 
% accuracy (CAA) > 80%, 50 bp < alignment length < 500 bp), along with the 
% 20 kinetic parameters (KP) used in the PCA analysis. Next, it selects the
% ones with the highest CAA for each pore ID with multiple potential 
% barcode hits. Finally, ...
%
% Case 1:
% ..., given a specific barcode ID (with the maximum frequency count), it 
% caluclates the PCA based on the 20 KP associated with this list of pores.
%
% Case 2: 
% ..., it randomly samles N pores from the appropriate barcode set, and 
% PCA-maps back as before. This way, the maixum number of N can be
% determined for "good" mapping. 
%
% This can be mapped back to PCA Fig. 6 in the main text for RPol
% verification. NOTE: Use MATLAB R2017a to run this code, older 
% version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function pca_mapback_pol_filter_matrix_norm(type)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> PCA map back filter matrix normalized start');
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
DATA = []; FILT = []; BAR = [];

% Generate 'experimental set' di rectory to hold the scatter plots.
cdir = 'pca_mapback';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Read filtered data matrix in.           
DATA = dlmread('data_matrix.txt');
            
% Analysis type 0 - all barcodes considered in PCA analysis.
if type == 0 

    % Retrieve the matrix of 20 kinetic properties only.
    FILT = DATA(:, 1:end);   

% Analysis type 1 - barcodes 1-32 considered in PCA analysis.
elseif type == 1

    % All barcodes in set 1.
    i1 = DATA(:,1)<33;
    FILT = DATA(i1, 1:end);

% Analysis type 2 - barcodes 33-64 considered in PCA analysis.
elseif type == 2

    % All barcodes in set 2.
    ii = DATA(:,1)>33;
    DATA2 = DATA(ii,:);
    i2 = DATA2(:,1)<65;
    FILT = DATA2(i2, 1:end);

% Analysis type 3 - barcodes 65-96 considered in PCA analysis.
elseif type == 3

    % All barcodes in set 3.
    i3 = DATA(:,1)>64;
    FILT = DATA(i3, 1:end);

else

    % Pick a specific barcode (maximum occurence, for example).
    ie = DATA(:,1) == type;
    FILT = DATA(ie, 1:end);

end

% Normalize counts by filtering out barcodes with at least O observations.
[N,edges] = histcounts(FILT(:,1), 96);
N

% Analysis type 0 - all barcodes considered in PCA analysis.
if type == 0 
    
    % Determine normalization parameters.
    Y = (sum(N(1:32))+sum(N(33:96)))/(length(N(1:32))+length(N(33:96)));
    
% Analysis type 1 - barcodes 1-32 considered in PCA analysis.
elseif type == 1

    Y = 0;
    
% Analysis type 2 - barcodes 33-64 considered in PCA analysis.
elseif type == 2
    
    Y = (sum(N(1:32))+sum(N(64:96)))/(length(N(1:32))+length(N(64:96)));
   
% Analysis type 3 - barcodes 65-96 considered in PCA analysis.
elseif type == 3
    
    Y = (sum(N(1:32))+sum(N(33:64)))/(length(N(1:32))+length(N(33:64)));
 
% Pick a specific barcode (maximum occurence, for example).    
else
    
    Y = 0;

end
    
% Determine average counts to deduct.
C = ceil(Y);
M = N-C;
M(M<0)=0;

% Build normalized barcode count histogram.
U=[];
for a = 1:length(M)
	for b = 1:M(a)
		U = [U; a];
	end
end

% Find unique barcode IDS.
[unique_bc_id, ia, ic] = unique(U);

% Construct nomalized kinetic property matrix.
for i = 1:length(unique_bc_id)
    
            in = FILT(:,1) == unique_bc_id(i);          
            BAR = vertcat(BAR, FILT(in, 4:end));
            
end
        
% Calculate principle components of data matrix.
opt = statset('pca'); opt.MaxIter = 1e+6;
[coeff, score, latent, tsquared, explained, mu] = pca(BAR, 'Algorithm', 'svd', 'Options', opt);
P = BAR*coeff;
PC1 = P(:,1); PC2 = P(:,2);

% Create 2D plot of PC1 vs. PC2.
scatter(PC1, PC2, 'o', 'LineWidth', 1);
grid;
title('Filtered PCA to Map Back');
xlabel('PC1');
ylabel('PC2');
axis([-4 8 -4 10]);
hold on;

BAR
size(BAR)
explained 

% Write data used in PCA to a matrix and save it.           
dlmwrite('pca_matrix.txt', BAR, 'delimiter', '\t', 'precision', '%.4f');

disp('--> PROCESSED DATA SET: data_matrix.txt'); 

% Save plot.
savefig([num2str(type), '_pca_mapback_pol_filter_norm.fig']);
print('-dbmp', [num2str(type), '_pca_mapback_pol_filter_norm.bmp']);

% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;   

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> PCA map back filter matrix normalized end');
fprintf('\n');

% End timer.
toc