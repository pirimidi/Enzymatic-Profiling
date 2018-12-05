%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: March 8, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters determined from a naive local 
% alignment (SW) algorithm. Then, for each 'pol-bar', it iterates through 
% the experimental structures and generates a consensus alignment for each 
% nanopore read, then determines a error profile for each ('pol-bar') 
% consensus. Finally, it generates whisker plots featuring consensus 
% alignment accuracy (%) versus number of full barcode iteration 
% ("align_copies") for each 'pol-bar' case displaying them on the same plot
% for comparison. NOTE: Use MATLAB R2017a to run this code, older version 
% might trigger errors. Vary 'pol', fix 'bar'.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function consensus_alignment

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Consensus alignment start');
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

% Set up regular expressions for inital template sequence.
expi = {'A', 'C', 'T', 'G'};
repi = {'N'; 'N'; 'N'; 'N'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DIRECTORY NAVIGATION SECTION');

% Navigate to 'pol-bar' MAT-file directory.
fas_dir = strcat(work_dir, '/data/mat_files');
cd(fas_dir);

disp(['--> IN DIRECTORY: ' fas_dir]);

% Read in FAS files for predefined barcodes (comp3, fv2, rep3).
bar_seqs = fastaread('barcode_list.fas');
                        
% Navigate to 'pol-bar' MAT-file directory.
mat_dir = strcat(work_dir, '/data/mat_files/cons');
cd(mat_dir);

disp(['--> IN DIRECTORY: ' mat_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('*pol6*_cons*.mat');
    
% Iterate through all 'pol' folders.
for e = 1:length(list) 
       
    % Create a data structure to store consensus object for each experiment. 
    cons_list = struct;
    
    % Load nanopore data structure per experimental set (pol-bar).
    load(list(e).name);
    disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
    % Iterate through all three barcodes.
    for bar = 1:length(barcodes)
        
        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cons_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cons_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cons_list(bar).experiment(exp).s;
                
                % For each nanopore read, calculate consensus sequence.
                for o = 1:length(pb)

                    % Define variables.
                    iteration = pb(o).iteration;
                    alignment = pb(o).alignment;
                    start = pb(o).start;
                    template = pb(o).lumped_template;

                    % Initiate array of fragment sequences.
                    seqs = {};
                    frag = {};
                    cseq = char.empty;

                    mc=1; fc=1;
                    fragment = regexprep(randseq(length(template(start(1):end))), expi, repi);

                    while ~strcmp(fragment, template(start(1):end))

                        cb = alignment(1+3*(mc-1));

                        if cb == 'A' || cb == 'C' || cb == 'T' || cb == 'G'
                            fragment(fc) = cb;
                            mc=mc+1; fc=fc+1;
                        else
                            mc=mc+1;
                        end
                    end

                    % First (iterative) fragment in alignment.
                    raw_1 = alignment(3+3*(0:mc-2));
                    raw_1l = regexprep(raw_1, '-', ''); % lump
                    seqs{1} = raw_1l;
                    frag{1} = raw_1;

                    % Do middle full (iterative) fragments; update counters.
                    if iteration > 2

                        for iter = 1:iteration-2

                            fc=1; sc=mc;
                            fragment = regexprep(randseq(length(template)), expi, repi);

                            while ~strcmp(fragment, template)

                                if(mc > length(alignment))
                                    break
                                end

                                cb = alignment(1+3*(mc-1));

                                if cb == 'A' || cb == 'C' || cb == 'T' || cb == 'G'
                                    fragment(fc) = cb;
                                    mc=mc+1; fc=fc+1;
                                else
                                    mc=mc+1;
                                end
                            end

                            % Middle (iterative) fragments in alignment.
                            raw_m = alignment(3+3*(sc-1:mc-2));
                            raw_ml = regexprep(raw_m, '-', ''); % lump
                            seqs{iter+1} = raw_ml;
                            frag{iter+1} = raw_m;

                        end
                    end

                    % Last (iterative) fragment in alignment.
                    raw_2 = alignment(3+3*(mc-1:length(alignment)-1));
                    raw_2l = regexprep(raw_2, '-', ''); % lump
                    seqs{iteration} = raw_2l;
                    frag{iteration} = raw_2;
                    
                    % Remove trailing empty fragment due to lumping.
                    if strcmp(seqs{end}, '') || strcmp(frag{end}, '')
                        seqs = seqs{1:end-1};
                        frag = frag{1:end-1};
                    end
                                        
                    % Correct for single fragment.
                    if iscell(frag)
                         
                        % Calculate consensus sequence and consensus identity.
                        if length(frag) == 1
                            cseq = seqs;
                            
                            % Update multialignment-based fields.
                            ma = struct;
                            profile = [];
                            wgtmatrix = {};
                            
                        elseif length(frag) == 2
                            [d_score, d_alignment, d_start] = swalign(seqs{1}, seqs{2});
                            [r, c] = size(d_alignment);
                                                        
                            for di = 1:c
                                if d_alignment(1+3*(di-1)) == '-'
                                    cseq(di) = d_alignment(3+3*(di-1));
                                else
                                    cseq(di) = d_alignment(1+3*(di-1));
                                end
                            end
                            
                            % Update multialignment-based fields.
                            ma = d_alignment;
                            profile = [];
                            wgtmatrix = {};                           

                        % Case for more than 3 iterations.
                        else
                            
                            ma = multialign(frag, 'ExistingGapAdjust', false, 'TerminalGapAdjust', true);
                            cseq = seqconsensus(ma, 'Alphabet', 'NT', 'Gaps', 'none');
                            [profile, symbols] = seqprofile(ma, 'Alphabet', 'NT', 'Counts', 'false', 'Gaps', 'none');
                            [wgtmatrix, handle] = seqlogo(ma, 'Displaylogo', 'false');
                      
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Calculate consensus identity to CORRECT template.
                        lumped_cseq = regexprep(cseq, expr, repl, 'ignorecase');
                        [c_score, c_alignment, c_start] = swalign(template, lumped_cseq);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Add sequence information to 'cons list'.
                        cons_list(bar).experiment(exp).s(o).lumped_fragments = seqs;
                        cons_list(bar).experiment(exp).s(o).alignment_fragments = frag;

                        cons_list(bar).experiment(exp).s(o).multialign = ma;
                        cons_list(bar).experiment(exp).s(o).consensus = cseq;
                        cons_list(bar).experiment(exp).s(o).lumped_consensus = lumped_cseq;
                        cons_list(bar).experiment(exp).s(o).profile = profile;
                        cons_list(bar).experiment(exp).s(o).weight_matrix = wgtmatrix;
                        
                        % Assign fields to CORRECT consensus alignment object.
                        cons_list(bar).experiment(exp).s(o).cons_score = c_score;
                        cons_list(bar).experiment(exp).s(o).cons_alignment = c_alignment;
                        cons_list(bar).experiment(exp).s(o).cons_alignment_length = length(c_alignment);
                        cons_list(bar).experiment(exp).s(o).cons_start = c_start;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Force alignment/consensus length to be at least 10 bp.
                        if length(c_alignment) >= 10
                        
                            % Initialize match, mismatch and deletion counters.
                            mat_count = 0;
                            mis_count = 0;  
                            ins_count = 0;
                            del_count = 0;

                            % Initialize match, mismatch and deletion index arrays.
                            mat_index = [];
                            mis_index = [];
                            ins_index = [];
                            del_index = [];  

                            % Compute number of matches, mismatchs and deletions.
                            for k = 1:length(c_alignment)

                                % An insertion.
                                %disp('--> Insertion');
                                if strcmp(c_alignment(1+3*(k-1)), '-')

                                    ins_count = ins_count + 1;
                                    ins_index = [ins_index k];

                                % A deletion.
                                %disp('--> Deletion');
                                elseif strcmp(c_alignment(3+3*(k-1)), '-')

                                    del_count = del_count + 1;
                                    del_index = [del_index k];

                                % A match or mismatch.    
                                else

                                    % Count number of matches.
                                    %disp('--> Match');
                                    if strcmp(c_alignment(1+3*(k-1)), c_alignment(3+3*(k-1)))

                                        mat_count = mat_count + 1;
                                        mat_index = [mat_index k];

                                    % Count number of mismatchs.
                                    %disp('--> mismatch');
                                    elseif not(strcmp(c_alignment(1+3*(k-1)), c_alignment(3+3*(k-1))))

                                        mis_count = mis_count + 1;
                                        mis_index = [mis_index k];

                                    end
                                end
                            end

                            % Calculate identity (%) of full alignment.
                            total = length(c_alignment(2, :));
                            match = length(find(c_alignment(2, :) == '|'));
                            identity = match/total*100; 

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Calculate consensus identity to two INCORRECT templates.                      
                            if strcmp(barcodes{bar}, 'comp3')                            
                                template_2 = bar_seqs(2).Sequence;
                                template_3 = bar_seqs(3).Sequence;

                            elseif strcmp(barcodes{bar}, 'fv2')
                                template_2 = bar_seqs(1).Sequence;
                                template_3 = bar_seqs(3).Sequence;

                            else
                                template_2 = bar_seqs(1).Sequence;
                                template_3 = bar_seqs(2).Sequence;

                            end

                            % Calculate identity (%) of full INCORRECT-2 alignment.
                            lumped_template_2 = regexprep(template_2, expr, repl, 'ignorecase');
                            [c_score_2, c_alignment_2, c_start_2] = swalign(lumped_template_2, cseq);

                            total_2 = length(c_alignment_2(2, :));
                            match_2 = length(find(c_alignment_2(2, :) == '|'));
                            identity_2 = match_2/total_2*100;

                            % Calculate identity (%) of full INCORRECT-3 alignment.
                            lumped_template_3 = regexprep(template_3, expr, repl, 'ignorecase');
                            [c_score_3, c_alignment_3, c_start_3] = swalign(lumped_template_3, cseq);

                            total_3 = length(c_alignment_3(2, :));
                            match_3 = length(find(c_alignment_3(2, :) == '|'));
                            identity_3 = match_3/total_3*100;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                            % Assign fields to INCORRECT-2 consensus alignment object.
                            cons_list(bar).experiment(exp).s(o).cons_score_2 = c_score_2;
                            cons_list(bar).experiment(exp).s(o).cons_alignment_2 = c_alignment_2;
                            cons_list(bar).experiment(exp).s(o).cons_alignment_length_2 = length(c_alignment_2);
                            cons_list(bar).experiment(exp).s(o).cons_start_2 = c_start_2;

                            % Assign fields to INCORRECT-3 consensus alignment object.
                            cons_list(bar).experiment(exp).s(o).cons_score_3 = c_score_3;
                            cons_list(bar).experiment(exp).s(o).cons_alignment_3 = c_alignment_3;
                            cons_list(bar).experiment(exp).s(o).cons_alignment_length_3 = length(c_alignment_3);
                            cons_list(bar).experiment(exp).s(o).cons_start_3 = c_start_3;

                            cons_list(bar).experiment(exp).s(o).cons_match_count = mat_count;
                            cons_list(bar).experiment(exp).s(o).cons_mismatch_count = mis_count;
                            cons_list(bar).experiment(exp).s(o).cons_insertion_count = ins_count;
                            cons_list(bar).experiment(exp).s(o).cons_deletion_count = del_count;

                            cons_list(bar).experiment(exp).s(o).cons_match_index = mat_index;
                            cons_list(bar).experiment(exp).s(o).cons_mismatch_index = mis_index;
                            cons_list(bar).experiment(exp).s(o).cons_insertion_index = ins_index; 
                            cons_list(bar).experiment(exp).s(o).cons_deletion_index = del_index;   

                            % Absolute index position in barcode.
                            if isempty(mat_index) 
                                cons_list(bar).experiment(exp).s(o).abs_cons_match_index = mat_index;
                            else
                                cons_list(bar).experiment(exp).s(o).abs_cons_match_index = mat_index + c_start(1);
                            end

                            if isempty(mis_index) 
                                cons_list(bar).experiment(exp).s(o).abs_cons_mismatch_index = mis_index;
                            else
                                cons_list(bar).experiment(exp).s(o).abs_cons_mismatch_index = mis_index + c_start(1);
                            end

                            if isempty(ins_index) 
                                cons_list(bar).experiment(exp).s(o).abs_cons_insertion_index = ins_index;
                            else
                                cons_list(bar).experiment(exp).s(o).abs_cons_insertion_index = ins_index + c_start(1);
                            end

                            if isempty(del_index) 
                                cons_list(bar).experiment(exp).s(o).abs_cons_deletion_index = del_index;
                            else
                                cons_list(bar).experiment(exp).s(o).abs_cons_deletion_index = del_index + c_start(1); 
                            end

                            cons_list(bar).experiment(exp).s(o).cons_identity = identity;
                            cons_list(bar).experiment(exp).s(o).cons_identity_2 = identity_2;
                            cons_list(bar).experiment(exp).s(o).cons_identity_3 = identity_3;

                            % INCORRECT template sequence information.
                            cons_list(bar).experiment(exp).s(o).raw_template_2 = template_2;
                            cons_list(bar).experiment(exp).s(o).raw_template_3 = template_3;      
                            cons_list(bar).experiment(exp).s(o).lumped_template_2 = lumped_template_2;
                            cons_list(bar).experiment(exp).s(o).lumped_template_3 = lumped_template_3;
                            
                        end
                    end                    
                end
            end
        end   
    end
    
    % Generate M-file containing all experimental statistics for each barcode.
    ss = strsplit(list(e).name, '.');
    pp = strsplit(ss{1}, '_');
    save([pp{1}, '_cons_align.mat'], 'cons_list');

end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Consensus alignment end');
fprintf('\n');

% End timer.
toc