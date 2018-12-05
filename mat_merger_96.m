%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: May 10, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a 'stats' and 'cons' MAT-file, one for 
% each experiment, containing all Genia-derived kinetic and consensus-
% derived parameters. Then, for each file, it iterates through the 
% experimental structures and merges 'stats' and 'cons' data into one file 
% based on cell ID. This way, the pore filtering algorithm can parse both 
% 'stats' and 'cons' derived data from one MAT structure. NOTE: Use MATLAB 
% R2017a to run this code, older version might trigger errors. 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function mat_merger_96

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> MAT merger 96 start');
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
mat_dir = strcat(work_dir, '/data/harvard/mat_files/');
cd(mat_dir);
disp(['--> IN DIRECTORY: ' mat_dir]);

% Read in a TXT file for predefined 96 naive barcodes.
bar_seqs = fastaread('unique_barcodes_50%.txt');

% Read in all 'experiment' folder names one-by-one.
cons = dir('*pol6*_cons_align*.mat');

% Read in all 'experiment' folder names one-by-one.
stats = dir('*pol6*_stats*.mat');
    
% Iterate through all 'cons' folders.
for e = 1:length(cons) 
       
    % Create a data structure to store consensus object for each experiment. 
    cons_list = struct; cells_list = struct; merge_list = struct;
    
    % Load nanopore data structure per experimental set (cons).
    cd(mat_dir);
    load(cons(e).name);
    disp(['--> MAT-FILE LOADED: ' cons(e).name]);
    load(stats(e).name);
    disp(['--> MAT-FILE LOADED: ' stats(e).name]);
        
    % Iterate through all 96 barcodes.
    for bar = 1:length(bar_seqs)
        
        % Iterate through all 'experiment' folders.
        for exp = 1:length(cons_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cons_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments and get cell ids. 
                ps = cells_list.experiment(exp).s;
                pb = cons_list(bar).experiment(exp).s;
                cell_id = extractfield(pb, 'cell_id');
                
                for ci = 1:length(cell_id)
                    in = find(strcmp({ps.cell_id}, cell_id{ci}) == 1);

                    % Merge to 'cons list' to 'cells_list'.
                    merge_list(bar).experiment(exp).s(ci) = ...
                    catstruct(cells_list.experiment(exp).s(in), ...
                    cons_list(bar).experiment(exp).s(ci));
                end
            end
        end   
    end
    
    % Generate M-file containing all experimental statistics for each barcode.
    ss = strsplit(cons(e).name, '.');
    pp = strsplit(ss{1}, '_');
    save([pp{1}, '_merge_96.mat'], 'merge_list', '-v7.3');

end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> MAT merger 96 end');
fprintf('\n');

% End timer.
toc