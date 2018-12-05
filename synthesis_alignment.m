%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: March 4, 2015.
%
% For: Quantification of de novo gene synthesis at the Church Lab - 
% Genetics Department, Harvard Medical School.
%
% Purpose: Given a DNA template and a list of *.fasta files, this program 
% iterates through all of them and extracts the base called sequences, 
% aligns them to the template and calculates: (1) length, (2) total count 
% and (3) error (insertion, deletion, mismatch) associated with each 
% single-pass alignment.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function synthesis_alignment

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic

fprintf('\n');
disp('--> Synthesis alignment start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working directory.
work_dir = pwd;

%-------------------------------------------------------------------------%
%                      MISEQ READS EXTRACTION SECTION                     %
%-------------------------------------------------------------------------%

disp('--> MISEQ READS EXTRACTION SECTION');

% Define current working directory.
cd 'extracted/60C/unfiltered';

% Read in AB, BC templates from FAS files.
AB_template = fastaread('AB.fas');
BC_template = fastaread('BC.fas');

% Read in all 'experiment' statistic text files one-by-one.
list = dir('*.fastq');     % all raw data
%list = dir('*.fasta');      % prefileterd data by Genious

% Read FASTA file into data structure.
alignment = {}; % clear contents (if any)
  
for i = 1:length(list)
  
    disp(['--> Processing file: ', list(i).name]);  
    
    %---------------------------------------------------------------------%
    %                     HANDLE RAW MISEQ READS                          %
    %---------------------------------------------------------------------%
    
    % Read FASTQ file into data structure.
    raw_data = {}; % clear contents (if any)
    raw_data = fastqread(list(i).name);    % all raw data
    %raw_data = fastaread(list(i).name);     % prefileterd data by Genious
    
    for j = 1:length(raw_data) 
        
        % Calculate template to single-pass, raw MiSeq read local alignment.
        [alignment(j).Score_AB, alignment(j).Alignment_AB] = swalign(AB_template.Sequence, raw_data(j).Sequence);
        [alignment(j).Score_BC, alignment(j).Alignment_BC] = swalign(BC_template.Sequence, raw_data(j).Sequence);
        
        % Check reverse complement as well.
        [alignment(j).Score_BA, alignment(j).Alignment_BA] = swalign(AB_template.Sequence, seqrcomplement(raw_data(j).Sequence));
        [alignment(j).Score_CB, alignment(j).Alignment_CB] = swalign(BC_template.Sequence, seqrcomplement(raw_data(j).Sequence));
        
        % Define constant arrays.
        score = [alignment(j).Score_AB alignment(j).Score_BC alignment(j).Score_BA alignment(j).Score_CB];
              
        % Find maximum alignment score and corresponsing sequence type.
        [M, I] = max(score);
        
        if I == 1
            alignment(j).Alignment = alignment(j).Alignment_AB;
            alignment(j).Type = 'AB';
        elseif I == 2
            alignment(j).Alignment = alignment(j).Alignment_BC;
            alignment(j).Type = 'BC';
        elseif I == 3
            alignment(j).Alignment = alignment(j).Alignment_BA;
            alignment(j).Type = 'AB';
        else
            alignment(j).Alignment = alignment(j).Alignment_CB;
            alignment(j).Type = 'BC';
        end            
        
        % Determine if it is a AB/BA or BC/CB sequence.
        alignment(j).Score = score(I);
        alignment(j).Length = length(raw_data(j).Sequence);
                            
        % Initialize match, mutation and deletion counters.
        mat_count = 0;
        mut_count = 0;  
        ins_count = 0;
        del_count = 0;
        
        % Initialize match, mutation and deletion index arrays.
        mat_index = [];
        mut_index = [];
        ins_index = [];
        del_index = [];  
        
        % Compute numer of matches, mutations and deletions.
        for k = 1:length(alignment(j).Alignment)

            % An insertion.
            %disp('--> Insertion');
            if strcmp(alignment(j).Alignment(1+3*(k-1)), '-')
                
                ins_count = ins_count + 1;
                ins_index = [ins_index k];
                
            % A deletion.
            %disp('--> Deletion');
            elseif strcmp(alignment(j).Alignment(3+3*(k-1)), '-')
                
                del_count = del_count + 1;
                del_index = [del_index k];
                
            % A match or mutation.    
            else
                                
                % Count number of matches.
                %disp('--> Match');
                if strcmp(alignment(j).Alignment(1+3*(k-1)), alignment(j).Alignment(3+3*(k-1)))

                    mat_count = mat_count + 1;
                    mat_index = [mat_index k];

                % Count number of mutations.
                %disp('--> Mutation');
                elseif not(strcmp(alignment(j).Alignment(1+3*(k-1)), alignment(j).Alignment(3+3*(k-1))))

                    mut_count = mut_count + 1;
                    mut_index = [mut_index k];

                end
            end
        end
        
        % Assign fields to alignment object.
        alignment(j).Match_count = mat_count;
        alignment(j).Mutation_count = mut_count;
        alignment(j).Insertion_count = ins_count;
        alignment(j).Deletion_count = del_count;
        
        alignment(j).Match_index = mat_index;
        alignment(j).Mutation_index = mut_index;
        alignment(j).Insertion_index = ins_index; 
        alignment(j).Deletion_index = del_index;    
           
    end    
    
    % Initialize counters.
    ab_count = 1;
    bc_count = 1;
    
    % Clear contents (if any).
    raw_AB_data = {};
    raw_BC_data = {};
    
    toc
    tic
    disp('--> Generate raw_NN_data sturctures');
    
    % Write the sequences to a FASTA-formatted file.
    for m = 1:length(raw_data) 
        
%         % Filter out reads with more than 10 mutations (junk).
%         if alignment(m).Mutation_count <= filter
       
            % Type 'AB' alignment.
            if strcmp(alignment(m).Type, 'AB')

                raw_AB_data(ab_count).Header = num2str(ab_count);
                raw_AB_data(ab_count).Sequence = raw_data(m).Sequence;

                raw_AB_data(ab_count).Score = alignment(m).Score;
                raw_AB_data(ab_count).Alignment = alignment(m).Alignment;
                raw_AB_data(ab_count).Length = alignment(m).Length;
                raw_AB_data(ab_count).Type = alignment(m).Type;

                raw_AB_data(ab_count).Match_count = alignment(m).Match_count;
                raw_AB_data(ab_count).Mutation_count = alignment(m).Mutation_count;
                raw_AB_data(ab_count).Insertion_count = alignment(m).Insertion_count;
                raw_AB_data(ab_count).Deletion_count = alignment(m).Deletion_count;

                raw_AB_data(ab_count).Match_index = alignment(m).Match_index;
                raw_AB_data(ab_count).Mutation_index = alignment(m).Mutation_index;
                raw_AB_data(ab_count).Insertion_index = alignment(m).Insertion_index;
                raw_AB_data(ab_count).Deletion_index = alignment(m).Deletion_index;

                % Increment counter.
                ab_count = ab_count + 1;

            % Type 'BC' alignment.    
            elseif strcmp(alignment(m).Type, 'BC')

                raw_BC_data(bc_count).Header = num2str(bc_count);
                raw_BC_data(bc_count).Sequence = raw_data(m).Sequence;

                raw_BC_data(bc_count).Score = alignment(m).Score;
                raw_BC_data(bc_count).Alignment = alignment(m).Alignment;
                raw_BC_data(bc_count).Length = alignment(m).Length;
                raw_BC_data(bc_count).Type = alignment(m).Type;

                raw_BC_data(bc_count).Match_count = alignment(m).Match_count;
                raw_BC_data(bc_count).Mutation_count = alignment(m).Mutation_count;
                raw_BC_data(bc_count).Insertion_count = alignment(m).Insertion_count;
                raw_BC_data(bc_count).Deletion_count = alignment(m).Deletion_count;

                raw_BC_data(bc_count).Match_index = alignment(m).Match_index;
                raw_BC_data(bc_count).Mutation_index = alignment(m).Mutation_index;
                raw_BC_data(bc_count).Insertion_index = alignment(m).Insertion_index;
                raw_BC_data(bc_count).Deletion_index = alignment(m).Deletion_index;

                % Increment counter.
                bc_count = bc_count + 1;

            end
%         end
    end
    
    %-------------------------------------------------------------------------%
    %                             FRAGMENT READS                              %
    %-------------------------------------------------------------------------%

    toc
    tic
    disp('--> Write 2 FASTA files');
    
    fastawrite([list(i).name, '_raw_AB_reads.txt'], raw_AB_data);
    fastawrite([list(i).name, '_raw_BC_reads.txt'], raw_BC_data);
    
    % Generate MAT-file structure for containing all these statistics.
    toc
    tic
    disp('--> Generate stats MAT file');

    m = matfile([list(i).name(1), '_stats.mat'], 'Writable', true);
    m.raw_AB_data = raw_AB_data; 
    m.raw_BC_data = raw_BC_data; 

    %-------------------------------------------------------------------------%
    %                            FIGURE GENERATION                            %
    %-------------------------------------------------------------------------%

    toc
    tic
    fprintf('\n');
    disp('--> Generate figures');

    % Create histogram of AB sequence length. 
    figure(1+(i-1)*10);
    histogram(extractfield(raw_AB_data, 'Length'), 1E2);   
    title('AB Sequence - Length Distribution');
    xlabel('Length (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_size_dist_AB.fig']);
    print('-dbmp', [list(i).name, '_size_dist_AB.bmp']); 
    disp('--> Length historgram created for AB sequence');

    % Display number of AB sequences.
    disp(['--> Total number of AB sequences: ', num2str(length(raw_AB_data))]);

    % Create histogram of AB sequence match per read. 
    figure(2+(i-1)*10);
    histogram(extractfield(raw_AB_data, 'Match_count'), 1E2);  
    title('AB Sequence - Match Per Read Distribution');
    xlabel('Match Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_match_dist_AB.fig']);
    print('-dbmp', [list(i).name, '_macth_dist_AB.bmp']); 
    disp('--> Match per read historgram created for AB sequence');

    % Create histogram of AB sequence mutation per read. 
    figure(3+(i-1)*10);
    histogram(extractfield(raw_AB_data, 'Mutation_count'), 1E1); 
    title('AB Sequence - Mutation Per Read Distribution');
    xlabel('Mutation Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_mutation_dist_AB.fig']);
    print('-dbmp', [list(i).name, '_mutation_dist_AB.bmp']); 
    disp('--> Mutation per read historgram created for AB sequence');

    % Create histogram of AB sequence insertion per read. 
    figure(4+(i-1)*10);
    histogram(extractfield(raw_AB_data, 'Insertion_count'), 2E1); 
    title('AB Sequence - Insertion Per Read Distribution');
    xlabel('Insertion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_insertion_dist_AB.fig']);
    print('-dbmp', [list(i).name, '_insertion_dist_AB.bmp']); 
    disp('--> Insertion per read historgram created for AB sequence');
    
    % Create histogram of AB sequence deletion per read. 
    figure(5+(i-1)*10);
    histogram(extractfield(raw_AB_data, 'Deletion_count'), 2E1); 
    title('AB Sequence - Deletion Per Read Distribution');
    xlabel('Deletion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_deletion_dist_AB.fig']);
    print('-dbmp', [list(i).name, '_deletion_dist_AB.bmp']); 
    disp('--> Deletion per read historgram created for AB sequence');

    % Create histogram of AB sequence match/mutation/insertion/deletion per position. 
    figure(6+(i-1)*10);
    ax1 = subplot(4,1,1); % top subplot
    ax2 = subplot(4,1,2); % mid-top subplot
    ax3 = subplot(4,1,3); % mid-bottom subplot
    ax4 = subplot(4,1,4); % bottom subplot

    title('AB Sequence - Position Changes');
    histogram(ax1, extractfield(raw_AB_data, 'Match_index'), length(AB_template.Sequence));
    xlabel(ax1, 'Match (position)');
    ylabel(ax1, 'Count (#)');
    xlim(ax1, [0 140]);
    histogram(ax2, extractfield(raw_AB_data, 'Mutation_index'), length(AB_template.Sequence));
    xlabel(ax2, 'Mutation (position)');
    ylabel(ax2, 'Count (#)');
    histogram(ax3, extractfield(raw_AB_data, 'Insertion_index'), length(AB_template.Sequence));
    xlabel(ax3, 'Insertion (position)');
    ylabel(ax3, 'Count (#)');
    histogram(ax4, extractfield(raw_AB_data, 'Deletion_index'), length(AB_template.Sequence));
    xlabel(ax4, 'Deletion (position)');
    ylabel(ax4, 'Count (#)');
    savefig([list(i).name, '_position_change_AB.fig']);
    print('-dbmp', [list(i).name, '_position_change_AB.bmp']); 
    disp('--> Position change historgram created for AB sequence');

    % Create histogram of BC sequence length. 
    figure(7+(i-1)*10);
    histogram(extractfield(raw_BC_data, 'Length'), 1E2);   
    title('BC Sequence - Length Distribution');
    xlabel('Length (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_size_dist_BC.fig']);
    print('-dbmp', [list(i).name, '_size_dist_BC.bmp']); 
    disp('--> Length historgram created for BC sequence');

    % Display number of BC sequences.
    disp(['--> Total number of BC sequences: ', num2str(length(raw_BC_data))]);

    % Create histogram of BC sequence match per read. 
    figure(8+(i-1)*10);
    histogram(extractfield(raw_BC_data, 'Match_count'), 1E2);  
    title('BC Sequence - Match Per Read Distribution');
    xlabel('Match Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_match_dist_BC.fig']);
    print('-dbmp', [list(i).name, '_macth_dist_BC.bmp']); 
    disp('--> Match per read historgram created for BC sequence');

    % Create histogram of BC sequence mutation per read. 
    figure(9+(i-1)*10);
    histogram(extractfield(raw_BC_data, 'Mutation_count'), 1E1); 
    title('BC Sequence - Mutation Per Read Distribution');
    xlabel('Mutation Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_mutation_dist_BC.fig']);
    print('-dbmp', [list(i).name, '_mutation_dist_BC.bmp']); 
    disp('--> Mutation per read historgram created for BC sequence');
    
    % Create histogram of BC sequence insertion per read. 
    figure(10+(i-1)*10);
    histogram(extractfield(raw_BC_data, 'Insertion_count'), 2E1); 
    title('BC Sequence - Insertion Per Read Distribution');
    xlabel('Insertion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_insertion_dist_BC.fig']);
    print('-dbmp', [list(i).name, '_insertion_dist_BC.bmp']); 
    disp('--> Insertion per read historgram created for BC sequence');
    
    % Create histogram of BC sequence deletion per read. 
    figure(11+(i-1)*10);
    histogram(extractfield(raw_BC_data, 'Deletion_count'), 2E1); 
    title('BC Sequence - Deletion Per Read Distribution');
    xlabel('Deletion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name, '_deletion_dist_BC.fig']);
    print('-dbmp', [list(i).name, '_deletion_dist_BC.bmp']); 
    disp('--> Deletion per read historgram created for BC sequence');

    % Create histogram of BC sequence match/mutation/deletion per position. 
    figure(12+(i-1)*10);
    ax1 = subplot(4,1,1); % top subplot
    ax2 = subplot(4,1,2); % mid-top subplot
    ax3 = subplot(4,1,3); % mid-bottom subplot
    ax4 = subplot(4,1,4); % bottom subplot

    title('BC Sequence - Position Changes');
    histogram(ax1, extractfield(raw_BC_data, 'Match_index'), length(BC_template.Sequence));
    xlabel(ax1, 'Match (position)');
    ylabel(ax1, 'Count (#)');
    xlim(ax1, [0 140]);
    histogram(ax2, extractfield(raw_BC_data, 'Mutation_index'), length(BC_template.Sequence));
    xlabel(ax2, 'Mutation (position)');
    ylabel(ax2, 'Count (#)');
    histogram(ax3, extractfield(raw_BC_data, 'Insertion_index'), length(BC_template.Sequence));
    xlabel(ax3, 'Insertion (position)');
    ylabel(ax3, 'Count (#)');
    histogram(ax4, extractfield(raw_BC_data, 'Deletion_index'), length(BC_template.Sequence));
    xlabel(ax4, 'Deletion (position)');
    ylabel(ax4, 'Count (#)');
    savefig([list(i).name, '_position_change_BC.fig']);
    print('-dbmp', [list(i).name, '_position_change_BC.bmp']); 
    disp('--> Position change historgram created for BC sequence');
    fprintf('\n');
    
    % Close all open figures.
    close all;

end

%-------------------------------------------------------------------------%
%                             OUTPUT ASSORTING                            %
%-------------------------------------------------------------------------%

disp('--> OUTPUT ASSORTING'); 
fprintf('\n');

% Move generated figures into 'stats_out' folder.
if ~exist('fig_out', 'dir')
  mkdir('fig_out');
end

movefile('*.fig', 'fig_out');
movefile('*.bmp', 'fig_out');

% Move generated FASTA-files into 'fas_out' folder.
if ~exist('fas_out', 'dir')
  mkdir('fas_out');
end

movefile('*reads.txt', 'fas_out');

% Move generated MAT-files into 'mat_out' folder.
if ~exist('mat_out', 'dir')
  mkdir('mat_out');
end

movefile('*.mat', 'mat_out');

% Navigate to working directory.
cd(work_dir);

% Stop timer.
toc

fprintf('\n');
disp('--> Synthesis alignment end');
fprintf('\n');
