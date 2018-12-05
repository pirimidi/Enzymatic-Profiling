%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: May 2, 2016.
%
% For: Quantification of de novo gene synthesis at the Church Lab - 
% Genetics Department, Harvard Medical School.
%
% Purpose: Given a data structure (*.mat) containing alignment information,
% this program generates statistics relevant to error types and prints the 
% output in a text file.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function int_stats

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic

fprintf('\n');
disp('--> Intermediate statistics start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working directory.
work_dir = pwd;

%-------------------------------------------------------------------------%
%                           STATISTICS SECTION                            %
%-------------------------------------------------------------------------%

disp('--> Statistics section');

% Navigate to *.mat alignment structure.
cd 'extracted/60C/unfiltered/mat_out';

% Read in all 'experiment' statistic text files one-by-one.
list = dir('*.mat');
  
% Open up text file for writing.
fid = fopen('stats_out.txt' ,'at');

% Initialize figure structure and counter.
f = {};
f_count = 0;
   
for h = 1:length(list)

    % Clear contents (if any).
    raw_AB_data = {};
    raw_BC_data = {};
    
    % Load alignment structures ['raw_N_data'] into workspace.
    disp(['--> Loading MAT-file: ', list(h).name]);
    load(list(h).name);

    % Raw data holder.
    rd = {raw_AB_data raw_BC_data};

    % Initialize counter structure.  
    c = {};

    for i = 1:length(rd)

        disp(['--> Processing structure: ', num2str(i)]);  
        
        % Increment figure counter.
        f_count = f_count + 1;
        
        % Initialize counters (for each fragment).
        bp_seq = 0;
        M = 0; D = 0; I = 0; 

        % Count total number of sequences. 
        tot_num = length(rd{i});

        % Count total number of bases sequenced. 
        bp_seq = sum(extractfield(rd{i}, 'Length'));

        % Count all errors (mutation + deletion + insertion).
        M = sum(extractfield(rd{i}, 'Mutation_count'));
        D = sum(extractfield(rd{i}, 'Deletion_count'));
        I = sum(extractfield(rd{i}, 'Insertion_count'));
        tot_err = M+I+D;

        % Count error frequency (per kB).
        err_freq = tot_err/bp_seq*1000;

        % Count perfect sequences.
        perfect_count = 0;
        perfect_index = [];

        for l = 1:length(rd{i})
            if rd{i}(l).Mutation_count == 0 && rd{i}(l).Deletion_count == 0 && rd{i}(l).Insertion_count == 0
                perfect_count = perfect_count + 1;
                perfect_index = [perfect_index l];
            end
        end

        disp('--> Count single deletions'); 

        % Count single deletions.
        single_dels = 0;
        multi_dels = 0;
        d_index = [];

        for q = 1:length(rd{i})

            % Pull out current deletion index array.
            d_index = rd{i}(q).Deletion_index;

            for r = 1:length(d_index)-1
                if d_index(r) + 1 ~= d_index(r+1)
                    single_dels = single_dels + 1;
                end
            end
        end

        % Calculate multiple deletions.
        multi_dels = D - single_dels;

        disp('--> Count single insertions');

        % Count single insertions.
        single_ins = 0;
        multi_ins = 0;
        i_index = [];

        for w = 1:length(rd{i})

            % Pull out current insertion index array.
            i_index = rd{i}(w).Insertion_index;

            for r = 1:length(i_index)-1
                if i_index(r) + 1 ~= i_index(r+1)
                    single_ins = single_ins + 1;
                end
            end
        end

        % Calculate multiple insertions.
        multi_ins = I - single_ins;            

        disp('--> Count type of mutations');

        % Count type of mutations.
        gc_count = 0; ga_count = 0; gt_count = 0;
        ag_count = 0; ac_count = 0; at_count = 0;
        ca_count = 0; ct_count = 0; cg_count = 0;
        ta_count = 0; tg_count = 0; tc_count = 0;

        m_index = [];

        for p = 1:length(rd{i})

            % Pull out current mutation index array.
            m_index = rd{i}(p).Mutation_index;

            for r = 1:length(m_index)

                % G>C mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'G') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'C')
                    gc_count = gc_count + 1;
                end

                % G>A mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'G') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'A')
                    ga_count = ga_count + 1;
                end

                % G>T mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'G') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'T')
                    gt_count = gt_count + 1;
                end

                % A>G mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'A') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'G')
                    ag_count = ag_count + 1;
                end

                % A>C mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'A') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'C')
                    ac_count = ac_count + 1;
                end

                % A>T mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'A') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'T')
                    at_count = at_count + 1;
                end

                % C>A mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'C') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'A')
                    ca_count = ca_count + 1;
                end

                % C>T mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'C') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'T')
                    ct_count = ct_count + 1;
                end

                % C>G mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'C') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'G')
                    cg_count = cg_count + 1;
                end

                % T>A mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'T') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'A')
                    ta_count = ta_count + 1;
                end

                % T>G mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'T') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'G')
                    tg_count = tg_count + 1;
                end

                % T>C mutation. 
                if strcmp(rd{i}(p).Alignment(1+3*(r-1)), 'T') && strcmp(rd{i}(p).Alignment(3+3*(r-1)), 'C')
                    tc_count = tc_count + 1;
                end

            end
        end    

        disp('--> Generate temporary count data structure');

        % Update current data structure with count values.  
        c(i).tot_num = tot_num;
        f(f_count).tot_num = tot_num;

        c(i).del = D;
        c(i).single_dels = single_dels;
        c(i).multi_dels = multi_dels;

        c(i).ins = I;
        c(i).single_ins = single_ins;
        c(i).multi_ins = multi_ins;

        c(i).mut = M;
        c(i).gc_count = gc_count; c(i).ga_count = ga_count; c(i).gt_count = gt_count;
        c(i).ag_count = ag_count; c(i).ac_count = ac_count; c(i).at_count = at_count;
        c(i).ca_count = ca_count; c(i).cg_count = cg_count; c(i).ct_count = ct_count;
        c(i).ta_count = ta_count; c(i).tg_count = tg_count; c(i).tc_count = tc_count;

        c(i).bp_seq = bp_seq; 
        c(i).tot_err = tot_err;
        c(i).err_freq = err_freq;
        c(i).per_seq = perfect_count;

    end

    % Print fragment statistics to file.
    fprintf(fid, '%s\n', list(h).name);
    struct2File(c, fid);
    fprintf(fid, '\n');
    
    %---------------------------------------------------------------------%
    %                      HISTOGRAM OF OCCURENCES                        %
    %---------------------------------------------------------------------%

    if(mod(h, 4) == 0)
        
        disp('--> Generate figures');

        % First set (h = 4).
        if h == 4
            
            % Total number per fragment.
            tot_num = [f(1).tot_num; f(2).tot_num; ...
                       f(3).tot_num; f(4).tot_num; ...
                       f(5).tot_num; f(6).tot_num; ...
                       f(7).tot_num; f(8).tot_num];
                   
        % Second set (h = 8).
        elseif h == 8
            
            % Total number per fragment.
            tot_num = [f(9).tot_num; f(10).tot_num; ...
                       f(11).tot_num; f(12).tot_num; ...
                       f(13).tot_num; f(14).tot_num; ...
                       f(15).tot_num; f(16).tot_num];
                   
        % Third set (h = 12).
        else
            
            % Total number per fragment.
            tot_num = [f(17).tot_num; f(18).tot_num; ...
                       f(19).tot_num; f(20).tot_num; ...
                       f(21).tot_num; f(22).tot_num; ...
                       f(23).tot_num; f(24).tot_num];
        end

        % Create normalized histogram of N sequence length. 
        figure(1);
        bar(tot_num);   
        title('Absolute Number of Sequences');
        xlabel('Cycles 1, 2, 3 & 5');
        ylabel('Absolute Count (#)');
        savefig([num2str(h/4), '_abs_dist.fig']);
        print('-dbmp', [num2str(h/4), '_abs_dist.bmp']); 
        disp('--> Absolute bargraph generated');

        figure(2);
        bar(tot_num/sum(tot_num));   
        title('Normalized Number of Sequences');
        xlabel('Cycles 1, 2, 3 & 5');
        ylabel('Normalized Count (#)');
        savefig([num2str(h/4), '_norm_dist.fig']);
        print('-dbmp', [num2str(h/4), '_norm_dist.bmp']); 
        disp('--> Normalized bargraph generated');

        % Close all open figures.
        close all;
    
    end
end
       
% Move generated figures into 'stats_out' folder.
if ~exist('bar_out', 'dir')
  mkdir('bar_out');
end

movefile('*.fig', 'bar_out');
movefile('*.bmp', 'bar_out');

% Close text file.
fclose(fid);

% Navigate to working directory.
cd(work_dir);

% Stop timer.
toc

fprintf('\n');
disp('--> Intermediate stats end');
fprintf('\n');
