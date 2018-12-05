%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: Aptil 5, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a 'stats' and 'cons' MAT-file, one for 
% each ('pol-bar') experiment, containing all Genia-derived kinetic and
% consensus-derived parameters. Then, for 'cons' each file, it iterates 
% through the experimental structures and merges 'stats' and 'cons' data
% into one file based on cell ID. This way, the pore filtering algorithm 
% can parse both 'stats' and 'cons' derived data from one MAT structure.
% NOTE: Use MATLAB R2017a to run this code, older version might trigger 
% errors. Vary 'pol', fix 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function mat_merger

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> MAT merger start');
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
                       
% Navigate to 'cons' MAT-file directory.
cons_dir = strcat(work_dir, '/data/mat_files/cons');
cd(cons_dir);
disp(['--> IN DIRECTORY: ' cons_dir]);

% Read in all 'experiment' folder names one-by-one.
cons = dir('*pol6*_cons*.mat');

% Navigate to 'stats' MAT-file directory.
stats_dir = strcat(work_dir, '/data/mat_files/stats');
cd(stats_dir);
disp(['--> IN DIRECTORY: ' stats_dir]);

% Read in all 'experiment' folder names one-by-one.
stats = dir('*pol6*_stats*.mat');
    
% Iterate through all 'pol' folders.
for e = 1:length(cons) 
       
    % Create a data structure to store consensus object for each experiment. 
    cons_list = struct; cells_list = struct; merge_list = struct;
    
    % Load nanopore data structure per experimental set (cons).
    cd(cons_dir);
    load(cons(e).name);
    disp(['--> MAT-FILE LOADED: ' cons(e).name]);
    cd(stats_dir);
    load(stats(e).name);
    disp(['--> MAT-FILE LOADED: ' stats(e).name]);
        
    % Iterate through all three barcodes.
    for bar = 1:length(barcodes)
        
        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cons_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cons_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and get cell ids. 
                ps = cells_list(bar).experiment(exp).s;
                pb = cons_list(bar).experiment(exp).s;
                cell_id = extractfield(pb, 'cell_id');
                
                for ci = 1:length(cell_id)
                    in = find(strcmp({ps.cell_id}, cell_id{ci}) == 1);

                    % Merge to 'cons list' to 'cells_list'.
                    merge_list(bar).experiment(exp).s(ci) = ...
                    catstruct(cells_list(bar).experiment(exp).s(in), ...
                    cons_list(bar).experiment(exp).s(ci));
                end
            end
        end   
    end
    
    % Generate M-file containing all experimental statistics for each barcode.
    ss = strsplit(cons(e).name, '.');
    pp = strsplit(ss{1}, '_');
    save([pp{1}, '_merge.mat'], 'merge_list');

end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> MAT merger end');
fprintf('\n');

% End timer.
toc