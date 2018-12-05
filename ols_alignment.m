%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: April 26, 2016.
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

function ols_alignment

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic

fprintf('\n');
disp('--> OLS alignment start');
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
cd 'test/DF0';

% Read in A, B, C, ..., I templates from FASTA files.
list = dir('M_*.fasta');

% Read in 'experiment' text file containing MiSeq raw reads.
raw_data = fastaread('Reads.fasta');
  
for i = 1:length(list)
  
    disp(['--> Processing file: ', list(i).name]);  
    
    %---------------------------------------------------------------------%
    %                     HANDLE RAW MISEQ READS                          %
    %---------------------------------------------------------------------%
    
    % Read FASTA file into data structure.
    alignment = {}; % clear contents (if any)
    template = fastaread(list(i).name);  
    
    for j = 1:length(raw_data) 
        
        % Calculate template to single-pass, raw MiSeq read local alignment.
        [alignment(i,j).Score, alignment(i,j).Alignment] = swalign(template.Sequence, raw_data(j).Sequence);
        
        % Check reverse complement as well.
        [alignment(i,j).Score_r, alignment(i,j).Alignment_r] = swalign(template.Sequence, seqrcomplement(raw_data(j).Sequence));
        
        % Define constant arrays.
        score = [alignment(i,j).Score alignment(i,j).Score_r];
              
        % Find maximum alignment score and corresponsing sequence type.
        [M, I] = max(score);
        
        if I == 1
            alignment(i,j).Alignment = alignment(i,j).Alignment;
            alignment(i,j).Type = 'n';
        else
            alignment(i,j).Alignment = alignment(i,j).Alignment_r;
            alignment(i,j).Type = 'r';
        end            
        
        % Determine score and length parameters.
        alignment(i,j).Score = score(I);
        alignment(i,j).Length = length(raw_data(j).Sequence);
                            
        % Initialize match, mutation, insertion and deletion counters.
        mat_count = 0;
        mut_count = 0;  
        ins_count = 0;
        del_count = 0;
        
        % Initialize match, mutation, insertion and deletion index arrays.
        mat_index = [];
        mut_index = [];
        ins_index = [];
        del_index = [];  
        
        % Compute numer of matches, mutations and deletions.
        for k = 1:length(alignment(i,j).Alignment)

            % An insertion.
            %disp('--> Insertion');
            if strcmp(alignment(i,j).Alignment(1+3*(k-1)), '-')
                
                ins_count = ins_count + 1;
                ins_index = [ins_index k];
                
            % A deletion.
            %disp('--> Deletion');
            elseif strcmp(alignment(i,j).Alignment(3+3*(k-1)), '-')
                
                del_count = del_count + 1;
                del_index = [del_index k];
                
            % A match or mutation.    
            else
                                
                % Count number of matches.
                %disp('--> Match');
                if strcmp(alignment(i,j).Alignment(1+3*(k-1)), alignment(i,j).Alignment(3+3*(k-1)))

                    mat_count = mat_count + 1;
                    mat_index = [mat_index k];

                % Count number of mutations.
                %disp('--> Mutation');
                elseif not(strcmp(alignment(i,j).Alignment(1+3*(k-1)), alignment(i,j).Alignment(3+3*(k-1))))

                    mut_count = mut_count + 1;
                    mut_index = [mut_index k];

                end
            end
        end
        
        % Assign fields to alignment object.
        alignment(i,j).Match_count = mat_count;
        alignment(i,j).Mutation_count = mut_count;
        alignment(i,j).Insertion_count = ins_count;
        alignment(i,j).Deletion_count = del_count;
        
        alignment(i,j).Match_index = mat_index;
        alignment(i,j).Mutation_index = mut_index;
        alignment(i,j).Insertion_index = ins_index; 
        alignment(i,j).Deletion_index = del_index;    
           
    end     
    
    % Initialize counter.
    count = 1;
    
    % Clear contents (if any).
    alignment_data = {};
    
    % Write the sequences to a FASTA-formatted file.
    for m = 1:length(raw_data) 
                   
        alignment_data(count).Header = num2str(count);
        alignment_data(count).Sequence = raw_data(m).Sequence;

        alignment_data(count).Score = alignment(i,m).Score;
        alignment_data(count).Alignment = alignment(i,m).Alignment;
        alignment_data(count).Length = alignment(i,m).Length;
        alignment_data(count).Type = alignment(i,m).Type;

        alignment_data(count).Match_count = alignment(i,m).Match_count;
        alignment_data(count).Mutation_count = alignment(i,m).Mutation_count;
        alignment_data(count).Insertion_count = alignment(i,m).Insertion_count;
        alignment_data(count).Deletion_count = alignment(i,m).Deletion_count;

        alignment_data(count).Match_index = alignment(i,m).Match_index;
        alignment_data(count).Mutation_index = alignment(i,m).Mutation_index;
        alignment_data(count).Insertion_index = alignment(i,m).Insertion_index;
        alignment_data(count).Deletion_index = alignment(i,m).Deletion_index;

        % Increment counter.
        count = count + 1;   
       
    end
    
    fastawrite([list(i).name(1:7), '_raw_reads.txt'], alignment_data);
    
    % Generate MAT-file structure for containing all these statistics.
    m = matfile([list(i).name(1:7), '_stats.mat'], 'Writable', true);
    m.alignment_data = alignment_data; 
       
    %-------------------------------------------------------------------------%
    %                            FIGURE GENERATION                            %
    %-------------------------------------------------------------------------%

    % Create histogram of fragment (A-I) length. 
    figure(1+(i-1)*10);
    histogram(extractfield(alignment_data, 'Length'), 1E2);   
    title([list(i).name(1:7), ' - Length Distribution']);
    xlabel('Length (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name(1:7), '_size_dist.fig']);
    print('-dbmp', [list(i).name(1:7), '_size_dist.bmp']); 
    disp(['--> Length historgram created for ', list(i).name(1:7), ' sequence']);

    % Display number of fragment (A-I).
    disp(['--> Total number of ', list(i).name(1:7), ' sequences: ', num2str(length(alignment_data))]);

    % Create histogram of fragment (A-I) match per read. 
    figure(2+(i-1)*10);
    histogram(extractfield(alignment_data, 'Match_count'), 1E2);  
    title([list(i).name(1:7), ' - Match Per Read Distribution']);
    xlabel('Match Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name(1:7), '_match_dist.fig']);
    print('-dbmp', [list(i).name(1:7), '_macth_dist.bmp']); 
    disp(['--> Match per read historgram created for ', list(i).name(1:7), ' sequence']);

    % Create histogram of fragment (A-I) mutation per read. 
    figure(3+(i-1)*10);
    histogram(extractfield(alignment_data, 'Mutation_count'), 1E2); 
    title([list(i).name(1:7), ' - Mutation Per Read Distribution']);
    xlabel('Mutation Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name(1:7), '_mutation_dist.fig']);
    print('-dbmp', [list(i).name, '_mutation_dist.bmp']); 
    disp(['--> Mutation per read historgram created for ', list(i).name(1:7), ' sequence']);

    % Create histogram of fragment (A-I) insertion per read. 
    figure(4+(i-1)*10);
    histogram(extractfield(alignment_data, 'Insertion_count'), 1E2); 
    title([list(i).name(1:7), ' - Insertion Per Read Distribution']);
    xlabel('Insertion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name(1:7), '_insertion_dist.fig']);
    print('-dbmp', [list(i).name(1:7), '_insertion_dist.bmp']); 
    disp(['--> Insertion per read historgram created for ', list(i).name(1:7), ' sequence']);
    
    % Create histogram of fragment (A-I) deletion per read. 
    figure(5+(i-1)*10);
    histogram(extractfield(alignment_data, 'Deletion_count'), 1E2); 
    title([list(i).name(1:7), ' - Deletion Per Read Distribution']);
    xlabel('Deletion Per Read (bp)');
    ylabel('Count (#)');
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
    grid;
    savefig([list(i).name(1:7), '_deletion_dist.fig']);
    print('-dbmp', [list(i).name(1:7), '_deletion_dist.bmp']); 
    disp(['--> Deletion per read historgram created for ', list(i).name(1:7), ' sequence']);

    % Create histogram of fragment (A-I) sequence match/mutation/insertion/deletion per position. 
    figure(6+(i-1)*10);
    ax1 = subplot(4,1,1); % top subplot
    ax2 = subplot(4,1,2); % mid-top subplot
    ax3 = subplot(4,1,3); % mid-bottom subplot
    ax4 = subplot(4,1,4); % bottom subplot

    title([list(i).name(1:7), ' - Position Changes']);
    histogram(ax1, extractfield(alignment_data, 'Match_index'), 1E2);
    xlabel(ax1, 'Match (position)');
    ylabel(ax1, 'Count (#)');
    histogram(ax2, extractfield(alignment_data, 'Mutation_index'), 1E2);
    xlabel(ax2, 'Mutation (position)');
    ylabel(ax2, 'Count (#)');
    histogram(ax3, extractfield(alignment_data, 'Insertion_index'), 1E2);
    xlabel(ax3, 'Insertion (position)');
    ylabel(ax3, 'Count (#)');
    histogram(ax4, extractfield(alignment_data, 'Deletion_index'), 1E2);
    xlabel(ax4, 'Deletion (position)');
    ylabel(ax4, 'Count (#)');
    savefig([list(i).name(1:7), '_position_change.fig']);
    print('-dbmp', [list(i).name(1:7), '_position_change.bmp']); 
    disp(['--> Position change historgram created for ', list(i).name(1:7), ' sequence']);
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
disp('--> OLS alignment end');
fprintf('\n');
