%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: May 9, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all alignment parameters generated by 
% 'consensus_alignment.m'. Then, for each 'pol-bar', it iterates through 
% the experimental structures and generates scatter plots featuring 
% consensus alignment accuracy (%) versus alignment length for each 
% 'pol-bar' case displaying them on the same plot for comparison. NOTE: 
% Use MATLAB R2017a to run this code, older version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function barcode_identity_96

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Barcode identity 96 start');
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

% Read in a TXT file for predefined 96 naive barcodes.
bar_seqs = fastaread('unique_barcodes_50%.txt');

% Read in all 'experiment' folder names one-by-one.
list = dir('*pol6*_cons_align*.mat');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> NORMALIZED SCATTER PLOTTING SECTION');

% Figure counter.
counter = 1;

% Generate 'experimental set' di rectory to hold the scatter plots.
cdir = 'barcode_identity_96';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Iterate through all MAT-files (currently only one).
 for e = 1:length(list) 

    % Load nanopore data structure per experimental set.
    load(list(e).name);
    disp(['--> MAT-FILE LOADED: ' list(e).name]);
    ss = strsplit(list(e).name, '.');

    % Iterate through all 'experiment' folders, here we have 8.
    for exp = 1:length(cons_list(1).experiment)

        % Define array container for fixed 'experiment' and grouping variable array.
        EX = []; GE = [];

        % Iterate through all 96 barcodes.
        for bar = 1:length(bar_seqs) 

            % Iterate through all experiments and plot desired field. 
            pb = cons_list(bar).experiment(exp).s;

            % Find filtered pores with valid 'cons_identity' value.
            i_1 = extractfield(pb, 'cons_identity') > 0;   
            pb_f1 = pb(i_1);
            
            % Only look at structures with elements (can be empty).
            if length(pb_f1) > 1
            
                % Find 'coi' values.
                norm = extractfield(pb_f1, 'cons_identity')./100;
                cell_id = extractfield(pb_f1, 'cell_id');

                % Generate container arrays for consensus identity, cell id and barcode id.
                EX = [EX; norm']; 
                GE = [GE; cell_id'];  
             
            end              

        end

        % Generate figure for each pol-bar combination.
        figure(counter);

        % Create dwell time (s) boxplot of all filtered events (exp).
        CategoricalScatterplot(EX, GE, 'Marker', 'o', 'FillMarker', true);
        grid;
        title('consensus identity vs. pore id');
        xlabel('pore id');
        ylabel('consensus identity');
        ylim([0 1.1]);
        xtickangle(45)
        legend(['pores = ', num2str(length(unique(GE)))], ['comparisons = ', num2str(length(EX))]);
        figureFullScreen(counter);

        % Save plot.
        fn = [num2str(exp), '_', 'bar_id_96'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;

        disp(['--> PROCESSED DATA SET: ' fn]); 
  
    end 
end    
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Barcode identity 96 end');
fprintf('\n');

% End timer.
toc