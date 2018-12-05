%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: May 10, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a MAT-file, containing all alignment 
% parameters generated by 'consensus_alignment_3.m'. Then, it iterates 
% through the experimental structures and collects all pore IDs, which are 
% categorized as a high-probability barcode hits (consensus alignment 
% accuracy (CAA) > 80%, 50 bp < alignment length < 500 bp). Finally, it 
% selects the ones with the highest CAA for each pore ID with multiple 
% potential barcode hits. NOTE: Use MATLAB R2017a to run this code, older 
% version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function barcode_all_unique_count_3(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Barcode all unique count 3 start');
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

% Read in a TXT file for predefined 3 naive barcodes.
bar_seqs = fastaread('barcode_list.fas');

% Read in all 'experiment' folder names one-by-one.
list = dir('*pol*_cons_align*.mat');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> NORMALIZED SCATTER PLOTTING SECTION');

% Figure counter.
counter = 1;

% Generate 'experimental set' di rectory to hold the scatter plots.
cdir = 'barcode_all_unique_count_3';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Define array container for all data.
AL = []; UC = 0;
    
% Iterate through all MAT-files (currently only one).
 for e = 1:length(list) 

    % Load nanopore data structure per experimental set.
    load(list(e).name);
    disp(['--> MAT-FILE LOADED: ' list(e).name]);
    ss = strsplit(list(e).name, '.');
        
    % Iterate through all 'experiment' folders, here we have 8.
    for exp = 1:length(cons_list(1).experiment)

        % Define array container for fixed 'experiment' and grouping variable array.
        CI = []; GE = []; BC = []; EX = [];
        ci = []; ge = {}; bc = []; ex = [];

        % Iterate through all 3 barcodes.
        for bar = 1:length(bar_seqs) 

            % Iterate through all experiments and plot desired field. 
            pb = cons_list(bar).experiment(exp).s;

            % Find filtered pores.
            i_1 = extractfield(pb, 'alignment_length') >= x_min;
            pb_f1 = pb(i_1);
            i_2 = extractfield(pb_f1, 'iteration') <= x_max;
            pb_f2 = pb_f1(i_2);
            i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
            pb_f3 = pb_f2(i_3);
            
            % Only look at structures with elements (can be empty).
            if length(pb_f3) > 1
            
                % Find 'coi' values.
                norm = extractfield(pb_f3, 'cons_identity')./100;
                cell_id = extractfield(pb_f3, 'cell_id');

                % Generate container arrays for consensus identity, cell id and barcode id.
                CI = [CI; norm']; 
                GE = [GE; cell_id']; 

                for b = 1:length(norm)
                    BC = [BC; bar'];
                end    
                
                for x = 1:length(norm)
                    EX = [EX; exp'];
                end                 
            end
        end

        % Filter out unique barcodes if more then one identified for a
        % single pore.
        
        [unique_cell_id, ia, ic] = unique(GE);
        
        for i = 1:length(unique_cell_id)
            in = find(strcmp(GE, unique_cell_id{i}) == 1);
            
            [ci(i), I] = max(CI(in)); % maximum consensus identity (%)
            rg = GE(in);
            ge{i} = rg{1}; % corresponding cell ID
            rb = BC(in);
            bc(i) = rb(I); % corresponding barcode ID
            re = EX(in);
            ex(i) = re(1); % corresponding experiment ID

        end
        
        % Update barcode array.
        UC = UC + length(unique_cell_id);
        AL = [AL, bc];
        
    end
    
    % Display data set (MAT file) processed.
    disp(['--> PROCESSED DATA SET: ' list(e).name]);
    
 end  

% Filter out unique barcodes if more then one identified for a
% single pore.
[unique_barcode_id, id, ie] = unique(AL);

% Create histogram for 'all' barcode hits (even multiple hits per pore).
histogram(AL, 3);
title('all unique barcode distribution');
xlabel('barcode id');
ylabel('count');
y_lim = ylim; axis([0 4 0 y_lim(2)+1]);
xticks(1:3);
xtickangle(45);
legend(['pores = ', num2str(UC), ...
    ' > barcodes = ', num2str(length(unique_barcode_id)), ...
    ' > total counts = ', num2str(length(AL))], ...
    'Location','northeastoutside');

% Save plot.
fn = 'bar_all_unique_count_3';
savefig([fn, '.fig']);
print('-dbmp', [fn, '.bmp']);
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Barcode all unique count 3 end');
fprintf('\n');

% End timer.
toc