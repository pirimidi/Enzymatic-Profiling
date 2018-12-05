%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: March 10, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all alignment parameters generated by 
% 'consensus_alignment.m'. Then, for each 'pol-bar', it iterates through 
% the experimental structures and creates 4 histograms, one for 
% match/mismatch/insertion/deletion per position in the consesnsus sequence
% displaying them on the same plot for comparison. NOTE: Use MATLAB R2017a
% to run this code, older version might trigger errors. Fix 'pol', vary 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function consensus_error_cons(type)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Consensus error pol start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
data_dir = strcat(work_dir, '\data\mat_files\cons');
cd(data_dir);

disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_cons*.mat');

% Figure counter.
counter = 1;

disp('--> ERROR HISTOGRAM PLOTTING SECTION');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         ERROR HISTOGRAM PLOTTING                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate 'experimental set' directory to hold the scatter plots.
cdir = ['consensus_error_cons_', type];

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Define array container for fixed 'bar' and grouping variable array.
MA = []; MU = []; IN = []; DE = [];
    
% Iterate through all 'pol' folders.
 for e = 1:length(list) 

    % Iterate through all three barcodes.
    for bar = 1:length(barcodes)  

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        ss = strsplit(list(e).name, '.');

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cons_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cons_list(bar).experiment(exp).s;

            % Find 'coi' values.
            ma = extractfield(pb, 'abs_cons_match_index');
            mu = extractfield(pb, 'abs_cons_mismatch_index');
            in = extractfield(pb, 'abs_cons_insertion_index');
            de = extractfield(pb, 'abs_cons_deletion_index');
            
            MA = [MA; ma']; MU = [MU; mu']; IN = [IN; in']; DE = [DE; de'];

        end
    end
 end
    
% Accunmulate all values to get maximum.
MM = [MA; MU; IN; DE];
        
% Generate figure for each pol-bar combination.
figure(counter);

% Create histogram of fragment (A-I) sequence match/mismatch/insertion/deletion per position. 
ax1 = subplot(4,1,1); % top subplot
ax2 = subplot(4,1,2); % mid-top subplot
ax3 = subplot(4,1,3); % mid-bottom subplot
ax4 = subplot(4,1,4); % bottom subplot

title('consensus error profile');

if type == 'p'

    histogram(ax1, MA, max(MA), 'Normalization', 'probability');
    ylabel(ax1, 'Probability');
    histogram(ax2, MU, max(MU), 'Normalization', 'probability');
    ylabel(ax2, 'Probability');
    histogram(ax3, IN, max(IN), 'Normalization', 'probability');
    ylabel(ax3, 'Probability');
    histogram(ax4, DE, max(DE), 'Normalization', 'probability');
    ylabel(ax4, 'Probability');

end

if type == 'c'

    histogram(ax1, MA, max(MA));
    ylabel(ax1, 'Count (#)');
    histogram(ax2, MU, max(MU));
    ylabel(ax2, 'Count (#)');
    histogram(ax3, IN, max(IN));
    ylabel(ax3, 'Count (#)');
    histogram(ax4, DE, max(DE));
    ylabel(ax4, 'Count (#)');
    
end

% Always the same for the plot.
xlabel(ax1, 'Match (position)');
xlim(ax1, [0 max(MM)]);
xlabel(ax2, 'Mismatch (position)');
xlim(ax2, [0 max(MM)]);
xlabel(ax3, 'Insertion (position)');
xlim(ax3, [0 max(MM)]);
xlabel(ax4, 'Deletion (position)');
xlim(ax4, [0 max(MM)]);

% Save plot.
fn = [type '_error_profile'];
savefig([fn, '.fig']);
print('-dbmp', [fn, '.bmp']);

% Update counter.
counter = counter + 1;

disp(['--> PROCESSED DATA SET: ' fn]);     
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Consensus error cons end');
fprintf('\n');

% End timer.
toc