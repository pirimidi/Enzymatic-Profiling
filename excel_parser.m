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

function excel_parser(folder_name)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> EXCEL parser start');
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

% Create a data structure to store cells object for each experiment. 
cells_list = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          DIRECTORY NAVIGATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
		for k = 1:r
		    CI{k} = RD{k,1};
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                                                                         %
		%                         CATEGORIES OF INTEREST                          %
		%                                                                         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		disp('--> CATEGORIES OF INTEREST SECTION');

		% Create an pre-populated cell array for keyset in hash table.
		coi = {'row (#)', 'col (#)', ...
			   'sequencing_end_time (seconds)', 'sequencing_lifetime (seconds)', 'lifetime_after_tag_flow (seconds)', 'len_starting_single_pore_run (#)', ...
			   'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
		       'align_copies (#)', ...
			   'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
			   'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)', ...
			   'level_call_A_counts (#)', 'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
			   'level_call_C_counts (#)', 'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
			   'level_call_G_counts (#)', 'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
			   'level_call_T_counts (#)', 'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)', ...
               'level_call_super_dwell_waiting_time_median (seconds)'};

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                                                                         %
		%                        DATA STRUCTURE BUILDING                          %
		%                                                                         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		disp('--> DATA STRUCTURE BUILDING SECTION');

		% Create a data structure to store numeric values for each 'coi' of "good" sequencing pore. 
		cells = struct;

		% Initialize cell counter.
		i = 1;

		% Iterate through all cell ids and build data structure with numeric values for each 'coi'. 
		for id_r = 2:length(CI)

			% Find value in 'cell annotations' corresponding to cell id/category pair and write it to data structure.
			cells(i).cell_id = CI{id_r};               
			cells(i).row = RD{id_r, find(strcmp(coi{1}, CA'))};
			cells(i).col = RD{id_r, find(strcmp(coi{2}, CA'))};

			% Sequencing time related parameters.
			cells(i).sequencing_end_time = RD{id_r, find(strcmp(coi{3}, CA'))};
			cells(i).sequencing_lifetime = RD{id_r, find(strcmp(coi{4}, CA'))};
			cells(i).lifetime_after_tag_flow = RD{id_r, find(strcmp(coi{5}, CA'))};
			cells(i).len_starting_single_pore_run = RD{id_r, find(strcmp(coi{6}, CA'))};

			% Ensamble-based base call and tag capture kinetics related parameters.
			cells(i).level_call_transition_rate = RD{id_r, find(strcmp(coi{7}, CA'))};
			cells(i).ttt_median = RD{id_r, find(strcmp(coi{8}, CA'))};
			cells(i).ttt_rate = RD{id_r, find(strcmp(coi{9}, CA'))};

			% Genia alignment related parameters.
			cells(i).align_copies = RD{id_r, find(strcmp(coi{10}, CA'))};
			
			% Heteropolymer (lumped) alignment.
			cells(i).align_align_length = RD{id_r, find(strcmp(coi{11}, CA'))};
			cells(i).align_num_delete = RD{id_r, find(strcmp(coi{12}, CA'))};
			cells(i).align_num_ident = RD{id_r, find(strcmp(coi{13}, CA'))};
			cells(i).align_num_insert = RD{id_r, find(strcmp(coi{14}, CA'))};
			cells(i).align_num_mismatch = RD{id_r, find(strcmp(coi{15}, CA'))};
			cells(i).align_procession_length = RD{id_r, find(strcmp(coi{16}, CA'))};
			cells(i).align_read_length = RD{id_r, find(strcmp(coi{17}, CA'))};

			% Homopolymer (unlumped) alignment.
			cells(i).align_homo_align_length = RD{id_r, find(strcmp(coi{18}, CA'))};
			cells(i).align_homo_num_delete = RD{id_r, find(strcmp(coi{19}, CA'))};
			cells(i).align_homo_num_ident = RD{id_r, find(strcmp(coi{20}, CA'))};
			cells(i).align_homo_num_insert = RD{id_r, find(strcmp(coi{21}, CA'))};
			cells(i).align_homo_num_mismatch = RD{id_r, find(strcmp(coi{22}, CA'))};
			cells(i).align_homo_procession_length = RD{id_r, find(strcmp(coi{23}, CA'))};
			cells(i).align_homo_read_length = RD{id_r, find(strcmp(coi{24}, CA'))};

			% Base-specific (A,C,G,T) kinetic parameters.
			cells(i).level_call_A_counts = RD{id_r, find(strcmp(coi{25}, CA'))};
			cells(i).level_call_A_k_cat_rate = RD{id_r, find(strcmp(coi{26}, CA'))};
			cells(i).level_call_A_k_off_rate = RD{id_r, find(strcmp(coi{27}, CA'))};
			cells(i).level_call_A_lower_bound = RD{id_r, find(strcmp(coi{28}, CA'))};
			cells(i).level_call_A_mean_dwell_time = RD{id_r, find(strcmp(coi{29}, CA'))};
			cells(i).level_call_A_ttt_median = RD{id_r, find(strcmp(coi{30}, CA'))};
			cells(i).level_call_A_ttt_rate = RD{id_r, find(strcmp(coi{31}, CA'))};
			cells(i).level_call_A_upper_bound = RD{id_r, find(strcmp(coi{32}, CA'))};

			cells(i).level_call_C_counts = RD{id_r, find(strcmp(coi{33}, CA'))};
			cells(i).level_call_C_k_cat_rate = RD{id_r, find(strcmp(coi{34}, CA'))};
			cells(i).level_call_C_k_off_rate = RD{id_r, find(strcmp(coi{35}, CA'))};
			cells(i).level_call_C_lower_bound = RD{id_r, find(strcmp(coi{36}, CA'))};
			cells(i).level_call_C_mean_dwell_time = RD{id_r, find(strcmp(coi{37}, CA'))};
			cells(i).level_call_C_ttt_median = RD{id_r, find(strcmp(coi{38}, CA'))};
			cells(i).level_call_C_ttt_rate = RD{id_r, find(strcmp(coi{39}, CA'))};
			cells(i).level_call_C_upper_bound = RD{id_r, find(strcmp(coi{40}, CA'))};

			cells(i).level_call_G_counts = RD{id_r, find(strcmp(coi{41}, CA'))};
			cells(i).level_call_G_k_cat_rate = RD{id_r, find(strcmp(coi{42}, CA'))};
			cells(i).level_call_G_k_off_rate = RD{id_r, find(strcmp(coi{43}, CA'))};
			cells(i).level_call_G_lower_bound = RD{id_r, find(strcmp(coi{44}, CA'))};
			cells(i).level_call_G_mean_dwell_time = RD{id_r, find(strcmp(coi{45}, CA'))};
			cells(i).level_call_G_ttt_median = RD{id_r, find(strcmp(coi{46}, CA'))};
			cells(i).level_call_G_ttt_rate = RD{id_r, find(strcmp(coi{47}, CA'))};
			cells(i).level_call_G_upper_bound = RD{id_r, find(strcmp(coi{48}, CA'))};

			cells(i).level_call_T_counts = RD{id_r, find(strcmp(coi{49}, CA'))};
			cells(i).level_call_T_k_cat_rate = RD{id_r, find(strcmp(coi{50}, CA'))};
			cells(i).level_call_T_k_off_rate = RD{id_r, find(strcmp(coi{51}, CA'))};
			cells(i).level_call_T_lower_bound = RD{id_r, find(strcmp(coi{52}, CA'))};
			cells(i).level_call_T_mean_dwell_time = RD{id_r, find(strcmp(coi{53}, CA'))};
			cells(i).level_call_T_ttt_median = RD{id_r, find(strcmp(coi{54}, CA'))};
			cells(i).level_call_T_ttt_rate = RD{id_r, find(strcmp(coi{55}, CA'))};
			cells(i).level_call_T_upper_bound = RD{id_r, find(strcmp(coi{56}, CA'))};
            
            cells(i).level_call_super_dwell_waiting_time_median = RD{id_r, find(strcmp(coi{57}, CA'))};

			% Increment counter.
			i = i + 1;
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                                                                         %
		%                          PARAMETER PLOTTING                             %
		%                                                                         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		disp('--> PARAMETER PLOTTING SECTION');

		% Define direcory to hold figures.
		if ~exist('plots', 'dir')
		  mkdir('plots');
        end
        
        % REMOVE after using it once!
        delete('pol6*');
        
		% Total number of cell ids.
		x = 1:length(CI)-1;

        for coi_id = 1:length(coi)
            
            % Ordinate name.
            coi_split = strsplit(coi{coi_id});
            coi_name = coi_split{1};
            ordinate_name = regexprep(coi_name, '_', '-');
            
            % Sort cells as desired for 'coi'.
            cells_sorted = nestedSortStruct2(cells, coi_name, -1);
            values = extractfield(cells_sorted, coi_name);
            keys = extractfield(cells_sorted, 'cell_id');

            % Sequencing end time vs. cell id plot.
            figure(coi_id);
            plot(x, values, 'rx', 'Color', col(col_id,:), 'MarkerSize', 2);
            grid;
            title([ordinate_name, ' vs cell id']);
            xlabel('cell id');
            ylabel(ordinate_name);
            legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
    		
            %xticks(x);
    		%xticklabels(keys);
    		%xtickangle(45);
            %figureFullScreen(coi_id, 'Full screen figure size');

            % Save plot.
            savefig([num2str(coi_id), '_', folder_name, '_', barcodes{bar}, '_', coi_name, '.fig']);
            print('-dbmp', [num2str(coi_id), '_', folder_name, '_', barcodes{bar}, '_', coi_name, '.bmp']);

            % Increment color ID.
            col_id = col_id + 1;
            
            % Reset color ID.
            if col_id == 21
                col_id = 1;
            end
        end
             
        % Read length histogram.
        figure(length(coi)+1);
        [ND, TD, RD] = xlsread('read_length_hist.csv');
        e=ND(:,1);
        x=ND(:,2);
        ed = [e; e(end)+1];
        histogram('BinEdges', ed','BinCounts', x);
        grid;
        title('Read Length Histogram');
        xlabel('Read Length (bp)');
		ylabel('Counts (#)');
        %axis([0 8 0 4*10^6]);
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+1, 'Full screen figure size');
        
		% Save plot.
		savefig([num2str(length(coi)+1), '_', folder_name, '_', barcodes{bar}, '_read_length_hist.fig']);
		print('-dbmp', [num2str(length(coi)+1), '_', folder_name, '_', barcodes{bar}, '_read_length_hist.bmp']);
        
        % Single pore lifetime histogram.
        figure(length(coi)+2);
        [ND, TD, RD] = xlsread('single_pore_hist.csv');
        e=ND(:,1);
        x=ND(:,2);
        y=ND(:,3);
        ed = [e; e(end)+1]*100; % Note, 1 rep = 100s.
        histogram('BinEdges',ed','BinCounts',x);
        grid;
        hold on;
        histogram('BinEdges',ed','BinCounts',y);
        set(gca,'YScale','log');
        title('Single Pore Lifetime Histogram');
        xlabel('Single pore lifetime (s)');
		ylabel('Counts (#)');
        %axis([0 8 0 4*10^6]);
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+2, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+2), '_', folder_name, '_', barcodes{bar}, '_single_pore_hist.fig']);
		print('-dbmp', [num2str(length(coi)+2), '_', folder_name, '_', barcodes{bar}, '_single_pore_hist.bmp']);

        % TTT histogram for each of the 4 tagged nucleotides.
        [ND, TD, RD] = xlsread('ttt_aggregate_hist.csv');
        el=ND(:,1);
        er=ND(:,2);
        ed = [el; er(end)];

        all=ND(:,3);
        c=ND(:,4);
        a=ND(:,5);
        t=ND(:,6);
        g=ND(:,7);

        % All tags.
        figure(length(coi)+3);
        histogram('BinEdges', ed', 'BinCounts', all);
        grid;
        set(gca,'YScale','log');
        title('Accumulative TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+3, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+3), '_', folder_name, '_', barcodes{bar}, '_ttt_all_hist.fig']);
		print('-dbmp', [num2str(length(coi)+3), '_', folder_name, '_', barcodes{bar}, '_ttt_all_hist.bmp']);

        % C-tag.
        figure(length(coi)+4);
        hc = histogram('BinEdges', ed', 'BinCounts', c);
        grid;
        y_lim = get(gca,'YLim');
        axis([0 8 0 y_lim(2)]);
        title('C-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        hc.FaceColor = 'blue';
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+4, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+4), '_', folder_name, '_', barcodes{bar}, '_ttt_c_hist.fig']);
		print('-dbmp', [num2str(length(coi)+4), '_', folder_name, '_', barcodes{bar}, '_ttt_c_hist.bmp']);

        % A-tag.
        figure(length(coi)+5);
        ha = histogram('BinEdges', ed', 'BinCounts', a);
        grid;
        axis([0 8 0 y_lim(2)]);
        title('A-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        ha.FaceColor = 'green';
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+5, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+5), '_', folder_name, '_', barcodes{bar}, '_ttt_a_hist.fig']);
		print('-dbmp', [num2str(length(coi)+5), '_', folder_name, '_', barcodes{bar}, '_ttt_a_hist.bmp']);

        % T-tag.
        figure(length(coi)+6);
        ht = histogram('BinEdges', ed', 'BinCounts', t);
        grid;
        axis([0 8 0 y_lim(2)]);
        title('T-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        ht.FaceColor = 'red';
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+6, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+6), '_', folder_name, '_', barcodes{bar}, '_ttt_t_hist.fig']);
		print('-dbmp', [num2str(length(coi)+6), '_', folder_name, '_', barcodes{bar}, '_ttt_t_hist.bmp']);

        % G-tag.
        figure(length(coi)+7);
        hg = histogram('BinEdges', ed', 'BinCounts', g);
        grid;
        axis([0 8 0 y_lim(2)]);
        title('G-tag TTT Histogram');
        xlabel('Time-to-Thread (ms)');
		ylabel('Counts (#)');
        hg.FaceColor = 'black';
        legend([folder_name, ' > ', barcodes{bar}, ' > ', regexprep(exp_dir, '_', '-'), ' (N = ', num2str(seq_pores), ')'], 'Location', 'northoutside');
        %figureFullScreen(length(coi)+7, 'Full screen figure size');
        
        % Save plot.
		savefig([num2str(length(coi)+7), '_', folder_name, '_', barcodes{bar}, '_ttt_g_hist.fig']);
		print('-dbmp', [num2str(length(coi)+7), '_', folder_name, '_', barcodes{bar}, '_ttt_g_hist.bmp']);
        
		% Move all figures to 'plots' directory.
		movefile('*.fig', 'plots');
		movefile('*.bmp', 'plots');
        
        % Close all open figures.
        close all;
                   
        % Add cells object to 'cells list'.
        cells_list(bar).experiment(exp).s = cells;
        
        % Move up a directory.
        cd('../');
    end   
end

% Move up a directory.
cd('../');
        
% Generate M-file containing all experimental statistics for each barcode.
save([folder_name, '_stats.mat'], 'cells_list');

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> EXCEL parser end');
fprintf('\n');

% End timer.
toc