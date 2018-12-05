 %--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: March 22, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program first randomly picks a MAT-file out of three 'pol'
% experimental sets, then randomly picks a 'pol-bar' experiment, followed
% a random pick of a pore with the correct barcode alignment. Then, it 
% remembers it's correct consensus identity (%) and alignment length. Also,
% it picks a random barcode from the unique barcode list (generated 
% previously with a 50% alignment identity to eachother). Subsequently,
% it aligns the "lumped" consesnus sequence to this unique and "lumped" 
% barcode using the Smith-Waterman algorithm. Finally, it plots the 
% randomly selected correct and incorrect consensus alignment % on the same 
% for comparison. Repeat this 10^3 times to build up statistics; and 
% finally determine X, Y cutoffs (e.g. ~50 bp; 60%). NOTE: Use MATLAB 
% R2017a to run this code, older version might trigger errors.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function lumped_consensus_scatter_collapsed_silico(num_iter)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Lumped consensus scatter collapsed silico start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};

% Set up regular expressions for "lumping" "stuttery" regions.
expr = {'A+', 'C+', 'T+', 'G+'};
repl = {'A'; 'C'; 'T'; 'G'};

% Create a data structure to store consensus object for each experiment. 
silico_list = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
fas_dir = strcat(work_dir, '\data\mat_files');
cd(fas_dir);

disp(['--> IN DIRECTORY: ' fas_dir]);

% Read in FAS files for predefined unique set of 96 barcodes with at most
% 50% overlap identity.
bar_seqs = fastaread('unique_barcodes_50%.txt');
                        
% Navigate to 'pol-bar' MAT-file directory.
mat_dir = strcat(work_dir, '/data/mat_files/cons');
cd(mat_dir);

disp(['--> IN DIRECTORY: ' mat_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('*pol6*_*cons*.mat');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCATTER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> NORMALIZED SCATTER PLOTTING SECTION');

% Generate 'experimental set' directory to hold the scatter plots.
cdir = 'lumped_consesnsus_scatter_collapsed_silico';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Define array container for 'coi' and grouping variable array.
COI = []; ge = []; 
COI_2 = []; ge_2 = []; 
       
% Initialize counter for alignments.
c = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate silico data n-times.
while (num_iter+1 ~= c)
    
    % Pick a random 'pol' folder.
    pol = randi(length(list)); 

    % Load nanopore data structure per experimental set (pol-bar).
    load(list(pol).name);
    %disp(['--> MAT-FILE LOADED: ' list(pol).name]);

    % Pick a random barcode.
    bar = randi(length(barcodes)); 

    % Pick a random experiment in this particular 'pol-bar' experimental set.
    exp = randi(length(cons_list(bar).experiment)); 

    % Only look at structures with elements (can be empty).
    if length(cons_list(bar).experiment(exp).s) > 1
               
        disp(['--> NUM ITER: ' num2str(c)]);
        
        % Pick a random pore in this 'pol-bar' experiment. 
        por = randi(length(cons_list(bar).experiment(exp).s));
        pb = cons_list(bar).experiment(exp).s(por);

        % Pull out CORRECT consesnsus sequence, identity and alignment length. 
        cseq = pb.consensus;
        norm = pb.cons_identity/100;
        align_length = pb.alignment_length;

        % Check for consensus length criteria.
        if(cseq > 10 & align_length > 50)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate random INCORRECT consensus sequence, identity and 
            % alignment length. 
            uni = randi(length(bar_seqs));
            template_2 = bar_seqs(uni).Sequence;

            % Calculate identity (%) of INCORRECT alignment.
            lumped_cseq = regexprep(cseq, expr, repl, 'ignorecase');
            lumped_template_2 = regexprep(template_2, expr, repl, 'ignorecase');
            [c_score_2, c_alignment_2, c_start_2] = swalign(lumped_template_2, lumped_cseq);

            total_2 = length(c_alignment_2(2, :));
            match_2 = length(find(c_alignment_2(2, :) == '|'));
            norm_2 = match_2/total_2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate container arrays for all experiments.
            COI = [COI; norm]; 
            COI_2 = [COI_2; norm_2];

            % Update concatinated grouping variable array.
            ge = [ge; 1];
            ge_2 = [ge_2; 2];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add consensus object to 'silico list'.
            silico_list(c).consensus = cseq;
            silico_list(c).lumped_consensus = lumped_cseq;
            silico_list(c).cons_identity = pb.identity;
            silico_list(c).align_length = align_length;
            
            silico_list(c).raw_template_2 = template_2;
            silico_list(c).lumped_template_2 = lumped_template_2;
            silico_list(c).cons_identity_2 = norm_2*100;  
            silico_list(c).cons_score_2 = c_score_2;
            silico_list(c).cons_alignment_2 = c_alignment_2;
            silico_list(c).cons_alignment_length_2 = length(c_alignment_2);
            silico_list(c).cons_start_2 = c_start_2;
            
            % Update counter.
            c = c + 1;
            
        end
    end 
end

% Generate figure for each pol-bar combination.
figure(1);

% Create dwell time (s) boxplot of all filtered events (exp).
% CategoricalScatterplot([COI; COI_2], [ge; ge_2], ...
%                        'WhiskerLine', false, ...
%                        'BoxColor', [0.8471 0.8627 0.8392], ...
%                        'WhiskerColor',[0.8235 0.7412 0.0392]);

CategoricalScatterplot([COI; COI_2], [ge; ge_2], ...
                       'Marker', 'o', ...
                       'MarkerSize', 2.5, ...
                       'FillMarker', false, ...
                       'spreadWidth', 0.1, ...
                       'boxWidth', 0.2, ...
                       'WhiskerLine', false, ...
                       'BoxColor', [0.8471 0.8627 0.8392], ...
                       'WhiskerColor',[0 0 1]);

% CategoricalScatterplot([COI; COI_2], [ge; ge_2], ...
%                        'Marker', 'o', ...
%                        'MarkerSize', 2, ...
%                        'FillMarker', true, ...
%                        'BoxColor', [0.8471 0.8627 0.8392], ...
%                        'BoxEdgeColor', [0 0 1], ...
%                        'MedianColor', [0 0 1], ...
%                        'WhiskerColor', [0.8235 0.7412 0.0392], ...
%                        'BoxAlpha', 0, ...
%                        'BoxLineStyle', '-', ...
%                        'WhiskerLine', false, ...
%                        'MedianLineStyle', '-', ...
%                        'WhiskerLineStyle', '-', ...
%                        'BoxLineWidth', 2, ...
%                        'MedianLineWidth', 2, ...
%                        'WhiskerLineWidth', 2);

grid;
title('in silico consensus identity vs. alignment type');
xlabel('alignment type');
ylabel('consensus identity');
ylim([0 1.1]);
legend('correct alignment', 'incorrect alignment');

% Save plot.
fn = 'lumped_consensus_scatter_collapsed_silico';
savefig([fn, '.fig']);
print('-dbmp', [fn, '.bmp']);
disp(['--> PROCESSED DATA SET: ' fn]);    
 
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Save silico alignment data to MAT file.
save([fn, '.mat'], 'silico_list');

% Navigate to working directory.
cd(work_dir);

disp('--> Lumped consensus scatter collapsed silico end');
fprintf('\n');

% End timer.
toc