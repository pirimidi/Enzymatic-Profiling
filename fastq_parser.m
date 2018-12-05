%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: February 6, 2017.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives an experimental folder name ('pol-barcode'),
% iterates through all sub-directories ('barcodes'), repeated experiment 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and collects categories of
% interests from the 'cell_annotations' EXCEL file. Then, using these para-
% meters, it plots user-defined correlations. NOTE: EXCEL file rows must be
% sorted by number of "good" pores at start ("Sort Largest to Smallest). 
% Use MATLAB R2017a to run this code, older version might trigger errors. 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function fastq_parser(folder_name)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> FASTQ parser start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
%barcodes = {'comp3', 'fv2', 'rep3'};
barcodes = {'fv2'};

% Color palette.
col = hsv(20); col_id = 1;

% Create a data structure to store consensus object for each experiment. 
cons_list = struct;

% Set up regular expressions for "lumping" "stuttery" regions.
expr = {'A+', 'C+', 'T+', 'G+'};
repl = {'A'; 'C'; 'T'; 'G'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to 'pol-bar' MAT-file directory.
fas_dir = strcat(work_dir, '\data\mat_files');
cd(fas_dir);

disp(['--> IN DIRECTORY: ' fas_dir]);

% Read in FAS files for predefined barcodes (comp3, fv2, rep3).
bar_seqs = fastaread('barcode_list.fas');

% Iterate through all three barcodes.
for bar = 1:length(barcodes) 

	disp('--> DIRECTORY NAVIGATION SECTION');

	% Navigate to 'pol-barcode' directory.
	data_dir = strcat(work_dir, '\data\a - individual_plots\', folder_name, '\', barcodes{bar});
	cd(data_dir);

	disp(['--> IN DIRECTORY: ' data_dir]);

	% Read in all 'experiment' folder names one-by-one.
	list = dir('*_HMS_*');

	% Iterate through all 'experiment' folders.
	for exp = 1:length(list) 

        % Clear data matrices.
        ND = {}; TD = {}; RD = {}; CA = {}; CI = {};
        
		% WARTORTLE generated experiment - Row 'BS' in 'cell_annotations.csv' file.
		if strcmp(list(exp).name(15), 'w')
			linenum=56;

		% PRIMEAPE generated experiment - Row 'S' in 'cell_annotations.csv' file.
		elseif strcmp(list(exp).name(15), 'p')
			linenum=58;

		else 
			disp('--> INCORRECT FOLDER PARSING!');
			exit;
        end
        
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                                                                         %
		%                       CELL ANNOTATIONS PARSING                          %
		%                                                                         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		disp('--> CELL ANNOTATIONS SECTION');

		% Navigate to 'experiment' directory.
		exp_dir = list(exp).name;
		cd(exp_dir);

		disp(['--> IN DIRECTORY: ' exp_dir]);

		% Get the number of functional sequencing pores from the 'metrics.txt' file. Pull out sorted EXCEL region 'A1:GAZ', where Z is 'seq_pores+1'.
		fid=fopen('metrics.txt');
		row = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', linenum-1);
		cow = row{1};
		n = regexp(cow, ['\d+\.?\d*'], 'match');
		seq_pores = str2num(cell2mat(n{1}(1)))
		fclose(fid);

		% Determine range of worksheet with pores of interest.
		xl_range = ['B1:GA' num2str(seq_pores+1)];

		% Read in data from cell annotations EXCEL file 'cell_annotations.csv'.
		[ND, TD, RD] = xlsread('cell_annotations.csv', xl_range);

		[r,c] = size(RD);

		% Generate category list from 'cell annotations'.
		for j = 1:c
		    CA{j} = RD{1,j};
		end

		% Generate cell id list from fist column of 'cell annotations'.
		for z = 1:r
		    CI{z} = RD{z,1};
        end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                                                                         %
		%                        DATA STRUCTURE BUILDING                          %
		%                                                                         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		disp('--> DATA STRUCTURE BUILDING SECTION');

		% Create a data structure to store numeric values for each 'coi' of "good" sequencing pore. 
		cons = struct;

        % Read in FASTQ text file containing heteropolymer raw reads.
        raw_data = fastqread('hetero_sequences.fastq');
  
		% Initialize cell counter.
		i = 1;
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                         %
        %               LUMPED ALIGNMENT AND CONSENSUS GENERATION                 %
        %                                                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('--> LUMPED ALIGNMENT AND CONSENSUS GENERATION SECTION');
            
		% Iterate through all cell ids and build data structure with numeric values for each 'coi'. 
		for id_r = 2:length(CI)
            
            % Find corresponding raw sequence from FASTQ file.
            index_c = strfind({raw_data.Header}, CI{id_r});
            index = find(not(cellfun('isempty', index_c)));
                       
            % Assign sequences.
            raw_template = bar_seqs(bar).Sequence;
            lumped_template = regexprep(raw_template, expr, repl, 'ignorecase'); % barcode
            lumped_sequence = raw_data(index).Sequence; % raw read
            lumped_template_con = [];            
            
            % Do consensus alignment as described in the algorithm 2.
            if length(lumped_sequence) > length(lumped_template)
                                
                % Calculate factor to concatinate barcode template 'f' times .
                iter = ceil(length(lumped_sequence) / length(lumped_template));
                
                % Generate container arrays for all experiments.
                for w = 1:iter              
                    lumped_template_con = [lumped_template_con lumped_template];                    
                end
                
                % Locally align template to single-pass, raw read using Smith-Waterman algorithm.
                [score, alignment, start] = swalign(lumped_template_con, lumped_sequence);  

                % Store values on data structure.
                cons(i).cell_id = CI{id_r};
                cons(i).header = raw_data(index).Header;
                cons(i).raw_sequence = raw_data(index).Sequence;
                cons(i).quality = raw_data(index).Quality;

                cons(i).raw_template = bar_seqs(bar).Sequence;
                cons(i).lumped_template = lumped_template;
                cons(i).lumped_sequence = lumped_sequence;
                cons(i).iteration = iter;
                cons(i).score = score;
                cons(i).alignment = alignment;
                cons(i).start = start;
                cons(i).alignment_length = length(alignment);

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
                for k = 1:length(alignment)

                    % An insertion.
                    %disp('--> Insertion');
                    if strcmp(alignment(1+3*(k-1)), '-')

                        ins_count = ins_count + 1;
                        ins_index = [ins_index k];

                    % A deletion.
                    %disp('--> Deletion');
                    elseif strcmp(alignment(3+3*(k-1)), '-')

                        del_count = del_count + 1;
                        del_index = [del_index k];

                    % A match or mismatch.    
                    else

                        % Count number of matches.
                        %disp('--> Match');
                        if strcmp(alignment(1+3*(k-1)), alignment(3+3*(k-1)))

                            mat_count = mat_count + 1;
                            mat_index = [mat_index k];

                        % Count number of mismatchs.
                        %disp('--> mismatch');
                        elseif not(strcmp(alignment(1+3*(k-1)), alignment(3+3*(k-1))))

                            mis_count = mis_count + 1;
                            mis_index = [mis_index k];

                        end
                    end
                end

                % Calculate identity (%) of full alignment.
                total = length(alignment(2, :));
                match = length(find(alignment(2, :) == '|'));
                identity = match/total*100; 

                % Assign fields to alignment object.
                cons(i).match_count = mat_count;
                cons(i).mismatch_count = mis_count;
                cons(i).insertion_count = ins_count;
                cons(i).deletion_count = del_count;

                % Relative index position in alignment.
                cons(i).match_index = mat_index;
                cons(i).mismatch_index = mis_index;
                cons(i).insertion_index = ins_index; 
                cons(i).deletion_index = del_index;   
                
                % Absolute index position in barcode.
                if isempty(mat_index) 
                    cons(i).abs_match_index = mat_index;
                else
                    cons(i).abs_match_index = mat_index + start(1);
                end
                
                if isempty(mis_index) 
                    cons(i).abs_mismatch_index = mis_index;
                else
                    cons(i).abs_mismatch_index = mis_index + start(1);
                end
                
                if isempty(ins_index) 
                    cons(i).abs_insertion_index = ins_index;
                else
                    cons(i).abs_insertion_index = ins_index + start(1);
                end
                
                if isempty(del_index) 
                    cons(i).abs_deletion_index = del_index;
                else
                    cons(i).abs_deletion_index = del_index + start(1); 
                end

                cons(i).identity = identity;

                % Increment counter.
                i = i + 1;
            end 

            % Add consensus object to 'cons list'.
            cons_list(bar).experiment(exp).s = cons;
            
        end
        
        % Move up a directory.
        cd('../');
        
    end   
end

% Move up a directory.
cd('../');
        
% Generate M-file containing all experimental statistics for each barcode.
save([folder_name, '_cons.mat'], 'cons_list');

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> FASTQ parser end');
fprintf('\n');

% End timer.
toc