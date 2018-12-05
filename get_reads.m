%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: June 3, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a list of *.fig files, this program iterates through all 
% of them and extracts the base called sequences from the tile names into a
% text file.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function get_reads

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

fprintf('\n');
disp('--> Get reads start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working directory.
work_dir = pwd;

%-------------------------------------------------------------------------%
%                          READS EXTRACTION SECTION                       %
%-------------------------------------------------------------------------%

disp('--> READS EXTRACTION SECTION');

% Define current working directory.
cd 'base_calls';

% Read in all 'experiment' statistic text files one-by-one.
list = dir('*.fig');

% Set up regular expressions for "lumping" "stuttery" regions.
expr = {'A+', 'C+', 'T+', 'G+'};
repl = {'A'; 'C'; 'T'; 'G'};
  
for i = 1:length(list)
  
    disp(['--> Processing file: ', list(i).name]);  

    % Get the title from opened MATLAB figure.
    open(list(i).name);
    h = get(gca, 'Title');
    t = get(h, 'String');
    
    index = findstr(cell2mat(t(1)), 'cell');    
    title = cell2mat(t(1));
    %cell = title(1:index+7);
    
    part = strsplit(list(i).name, '_');
    cell = [part{1} '_' part{2}]; 
       
    % Get the template read from the first MATLAB figure. No need to
    % iterate on this, since it is always the same per script run.
    if i == 1
        raw_template = cell2mat(regexprep(t(4), '[-]', ''));
        lumped_template = regexprep(raw_template, expr, repl, 'ignorecase'); 
    end
    
    %---------------------------------------------------------------------%
    %                       HANDLE RAW READS                              %
    %---------------------------------------------------------------------%
    
    % Remove '-' character from nanopore read and convert it to a string
    % (from cell) data type.
    raw_read = cell2mat(regexprep(t(2), '[-]', ''));   
    
    % Write this into FASTA-formatted structure keeping "stutter" as is.
    raw_data(i).Sequence = raw_read;
    raw_data(i).Header =  cell;
     
    % Calculate template to single-pass, raw nanopore read global/local 
    % alignment.
    raw_data(i).Global = nwalign(raw_template, raw_read);
    raw_data(i).Local = swalign(raw_template, raw_read);
        
    %---------------------------------------------------------------------%
    %                      HANDLE LUMPED READS                            %
    %---------------------------------------------------------------------%
    
    % Collapse homopolymer runs into single nucleotides.
    lumped_read = regexprep(raw_read, expr, repl, 'ignorecase');   
        
    % Write this into FASTA-formatted structure keeping "stutter" as is.
    lumped_data(i).Sequence = lumped_read;
    lumped_data(i).Header =  cell;
    
    % Calculate template to single-pass, lumped nanopore read global/local 
    % alignment.
    lumped_data(i).Global = nwalign(lumped_template, lumped_read);
    lumped_data(i).Local = swalign(lumped_template, lumped_read);
        
    % Close opened MATLAB figure.
    close;
end

% Write the sequences to a FASTA-formatted file.
fastawrite('../raw_reads.txt', raw_data);
fastawrite('../lumped_reads.txt', lumped_data);

%-------------------------------------------------------------------------%
%                     SINGLE-PASS SEQUENCE ALIGNMENT                      %
%-------------------------------------------------------------------------%

% Sort raw reads by descending global/local alignment score.
raw_global_sorted = nestedSortStruct(raw_data, 'Global');
raw_local_sorted = nestedSortStruct(raw_data, 'Local');

% Sort lumped reads by descending global/local alignment score.
lumped_global_sorted = nestedSortStruct(lumped_data, 'Global');
lumped_local_sorted = nestedSortStruct(lumped_data, 'Local');

% Calculate template to single-pass nanopore read global/local alignment.
[raw_global_score_i, raw_global_alignment_i] = nwalign(raw_template,  raw_global_sorted(length(list)).Sequence, 'Showscore', 'true');
savefig('../raw_global_path_i.fig');
[lumped_global_score_i, lumped_global_alignment_i] = nwalign(lumped_template, lumped_global_sorted(length(list)).Sequence, 'Showscore', 'true');
savefig('../lumped_global_path_i.fig');

[raw_local_score_i, raw_local_alignment_i] = swalign(raw_template, raw_local_sorted(length(list)).Sequence, 'Showscore', 'true');
savefig('../raw_local_path_i.fig');
[lumped_local_score_i, lumped_local_alignment_i] = swalign(lumped_template, lumped_local_sorted(length(list)).Sequence, 'Showscore', 'true');
savefig('../lumped_local_path_i.fig');

%-------------------------------------------------------------------------%
%                             TERMINAL DISPLAY                            %
%-------------------------------------------------------------------------%

% Display single-pass sequences for both "raw" and "lumped" global/local 
% alignments.
fprintf('\n\n');
disp('--> BEST GLOBAL SINGLE-PASS ALIGNMENT');

raw_template
disp(['--> Best raw score: ' num2str(raw_global_score_i)]);
total = length(raw_global_alignment_i(2, :));
match = length(find(raw_global_alignment_i(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
disp(['--> Cell ID: ' raw_global_sorted(length(list)).Header]);
raw_global_alignment_i
showalignment(raw_global_alignment_i);

lumped_template
disp(['--> Best lumped score: ' num2str(lumped_global_score_i)]);
total = length(lumped_global_alignment_i(2, :));
match = length(find(lumped_global_alignment_i(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
disp(['--> Cell ID: ' lumped_global_sorted(length(list)).Header]);
lumped_global_alignment_i
showalignment(lumped_global_alignment_i);
fprintf('\n');

disp('--> BEST LOCAL SINGLE-PASS ALIGNMENT');

raw_template
disp(['--> Best raw score: ' num2str(raw_local_score_i)]);
total = length(raw_local_alignment_i(2, :));
match = length(find(raw_local_alignment_i(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
disp(['--> Cell ID: ' raw_local_sorted(length(list)).Header]);
raw_local_alignment_i
showalignment(raw_local_alignment_i);

lumped_template
disp(['--> Best lumped score: ' num2str(lumped_local_score_i)]);
total = length(lumped_local_alignment_i(2, :));
match = length(find(lumped_local_alignment_i(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
disp(['--> Cell ID: ' lumped_local_sorted(length(list)).Header]);
lumped_local_alignment_i
showalignment(lumped_local_alignment_i);

%-------------------------------------------------------------------------%
%                        SEQUENCE ALIGNMENT OUTPUT                        %
%-------------------------------------------------------------------------%

for j = 1:length(list)    
    if j == 1
        fprintf('\n');
        disp('--> RAW GLOBAL SINGLE-PASS ALIGNMENT');
        raw_template          
    end

    % Calculate current template to single-pass nanopore read alignment.
    [raw_global_score_c, raw_global_alignment_c] = nwalign(raw_template,  raw_global_sorted(length(list)-(j-1)).Sequence);   
    disp(['--> Raw score: ' num2str(raw_global_score_c)]);
    total = length(raw_global_alignment_c(2, :));
    match = length(find(raw_global_alignment_c(2, :) == '|'));
    disp(['--> Identity: ' num2str(match/total*100) '%']);
    disp(['--> Cell ID: ' raw_global_sorted(length(list)-(j-1)).Header]);
    raw_global_alignment_c    
end

for k = 1:length(list)    
    if k == 1
        fprintf('\n');
        disp('--> LUMPED GLOBAL SINGLE-PASS ALIGNMENT');
        lumped_template    
    end

    % Calculate current template to single-pass nanopore read alignment.
    [lumped_global_score_c, lumped_global_alignment_c] = nwalign(lumped_template, lumped_global_sorted(length(list)-(k-1)).Sequence);
    disp(['--> Lumped score: ' num2str(lumped_global_score_c)]);
    total = length(lumped_global_alignment_c(2, :));
    match = length(find(lumped_global_alignment_c(2, :) == '|'));
    disp(['--> Identity: ' num2str(match/total*100) '%']);
    disp(['--> Cell ID: ' lumped_global_sorted(length(list)-(k-1)).Header]);
    lumped_global_alignment_c  
end

for l = 1:length(list)    
    if l == 1   
        fprintf('\n');
        disp('--> RAW LOCAL SINGLE-PASS ALIGNMENT');
        raw_template         
    end

    % Calculate current template to single-pass nanopore read alignment.
    [raw_local_score_c, raw_local_alignment_c] = swalign(raw_template, raw_local_sorted(length(list)-(l-1)).Sequence);
    disp(['--> Raw score: ' num2str(raw_local_score_c)]);
    total = length(raw_local_alignment_c(2, :));
    match = length(find(raw_local_alignment_c(2, :) == '|'));
    disp(['--> Identity: ' num2str(match/total*100) '%']); 
    disp(['--> Cell ID: ' raw_local_sorted(length(list)-(l-1)).Header]);
    raw_local_alignment_c
 
end

for m = 1:length(list)    
    if m == 1  
        fprintf('\n');
        disp('--> LUMPED LOCAL SINGLE-PASS ALIGNMENT');
        lumped_template       
    end

    % Calculate current template to single-pass nanopore read alignment.
    [lumped_local_score_c, lumped_local_alignment_c] = swalign(lumped_template, lumped_local_sorted(length(list)-(m-1)).Sequence);
    disp(['--> Lumped score: ' num2str(lumped_local_score_c)]);
    total = length(lumped_local_alignment_c(2, :));
    match = length(find(lumped_local_alignment_c(2, :) == '|'));
    disp(['--> Identity: ' num2str(match/total*100) '%']);
    disp(['--> Cell ID: ' lumped_local_sorted(length(list)-(m-1)).Header]);    
    lumped_local_alignment_c    
end

%-------------------------------------------------------------------------%
%                       MULTIPLE SEQUENCE ALIGNMENT                       %
%-------------------------------------------------------------------------%

% Align multiple sequences using progressive method.
raw_reads = fastaread('../raw_reads.txt');
mar = multialign(raw_reads);
seqalignviewer(mar, 'SeqHeaders', raw_reads);

lumped_reads = fastaread('../lumped_reads.txt');
mal = multialign(lumped_reads);
seqalignviewer(mal, 'SeqHeaders', lumped_reads);

% Finally calculate consensus sequence based on this multiply aligned
% nucleotide sequences.
seqlogo(mar);
[raw_consensus, raw_score] = seqconsensus(mar);

seqlogo(mal);
[lumped_consensus, lumped_score] = seqconsensus(mal);

% Calculate template to consesnsus nanopore read global/local alignment.
[raw_global_score, raw_global_alignment] = nwalign(raw_template, raw_consensus, 'Showscore', 'true');
savefig('../raw_global_path.fig');
[lumped_global_score, lumped_global_alignment] = nwalign(lumped_template, lumped_consensus, 'Showscore', 'true');
savefig('../lumped_global_path.fig');

[raw_local_score, raw_local_alignment] = swalign(raw_template, raw_consensus, 'Showscore', 'true');
savefig('../raw_local_path.fig');
[lumped_local_score, lumped_local_alignment] = swalign(lumped_template, lumped_consensus, 'Showscore', 'true');
savefig('../lumped_local_path.fig');

%-------------------------------------------------------------------------%
%                             TERMINAL DISPLAY                            %
%-------------------------------------------------------------------------%

% Display consensus sequences for both "raw" and "lumped" global/local 
% alignments.
fprintf('\n');
disp('--> GLOBAL CONSENSUS ALIGNMENT');

raw_template
disp(['--> Raw score: ' num2str(raw_global_score)]);
total = length(raw_global_alignment(2, :));
match = length(find(raw_global_alignment(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
raw_global_alignment
showalignment(raw_global_alignment);

lumped_template
disp(['--> Lumped score: ' num2str(lumped_global_score)]);
total = length(lumped_global_alignment(2, :));
match = length(find(lumped_global_alignment(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
lumped_global_alignment
showalignment(lumped_global_alignment);

fprintf('\n');
disp('--> LOCAL CONSENSUS ALIGNMENT');

raw_template
disp(['--> Raw score: ' num2str(raw_local_score)]);
total = length(raw_local_alignment(2, :));
match = length(find(raw_local_alignment(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
raw_local_alignment
showalignment(raw_local_alignment);

lumped_template
disp(['--> Lumped score: ' num2str(lumped_local_score)]);
total = length(lumped_local_alignment(2, :));
match = length(find(lumped_local_alignment(2, :) == '|'));
disp(['--> Identity: ' num2str(match/total*100) '%']);
lumped_local_alignment
showalignment(lumped_local_alignment);

%-------------------------------------------------------------------------%
%                             OUTPUT ASSORTING                            %
%-------------------------------------------------------------------------%

% Navigate to working directory.
cd(work_dir);

% Move statistics files into 'stats_out' folder.
if ~exist('stats_out', 'dir')
  mkdir('stats_out');
end

movefile('*.txt', 'stats_out');

% Move generated figures into 'stats_out' folder.
if ~exist('figus_out', 'dir')
  mkdir('figus_out');
end

movefile('*.fig', 'figus_out');

disp('--> Get reads end');
fprintf('\n');
