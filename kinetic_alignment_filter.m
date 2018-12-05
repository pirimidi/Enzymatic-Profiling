%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: April 7, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters parsed from the 'cell_annotations' 
% EXCEL file. Then, for each 'pol-bar', it iterates through the experimental 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates normalized 
% whisker plots featuring alignment accuracy (%) versus number of full 
% barcode iteration ("align_copies") for each filtered 'pol-bar' case 
% displaying them on the same plot for comparison. NOTE: Use MATLAB R2017a 
% to run this code, older version might trigger errors. Vary 'pol', fix 'bar'.
%
% DEFINITIONS:      
%
% (1) heteropolymer alignment accuracy (%) = align_num_ident / align_align_length
% (2) homopolymer alignment accuracy (%) = align_homo_num_ident / align_homo_align_length
% 
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_alignment_filter(x_min, x_max, y_min)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic alignment filter start');
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

% Navigate to 'pol-bar' MAT-file directory.
data_dir = strcat(work_dir, '\data\mat_files\merge');
cd(data_dir);

disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_merge*.mat');

% Figure counter.
counter = 1;

disp('--> NORMALIZED WHISKER PLOTTING SECTION');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              WHISKER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate 'experimental set' directory to hold the whisker plots.
cdir = 'kinetic_alignment_filter';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Iterate through all three barcodes.
for bar = 1:length(barcodes)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for heteropolymer deletions.
    coi_name = 'align_num_delete';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                    
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;
                
                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
            
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for heteropolymer identicals.
    coi_name = 'align_num_ident';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
            
            end          
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for heteropolymer insertions.
    coi_name = 'align_num_insert';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

               % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
            
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for heteropolymer mismatches.
    coi_name = 'align_num_mismatch';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
                
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for homopolymer deletions.
    coi_name = 'align_homo_num_delete';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_homo_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
                
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for homopolymer identicals.
    coi_name = 'align_homo_num_ident';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_homo_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
                
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filter ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for homopolymer insertions.
    coi_name = 'align_homo_num_insert';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_homo_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
                
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filter ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for homopolymer mismatches.
    coi_name = 'align_homo_num_mismatch';
    ordinate_name = regexprep(coi_name, '_', '-');

    % Generate figure for each pol-bar combination.
    figure(counter);
    
    % Define array container for fixed 'bar'.
    BAR = [];

    % Iterate through all 'pol' folders.
    for e = 1:length(list) 

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        
        % Define array container for 'coi'.
        COI = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(merge_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(merge_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = merge_list(bar).experiment(exp).s;

                % Find filtered pores.
                i_1 = extractfield(pb, 'alignment_length') >= x_min;
                pb_f1 = pb(i_1);
                i_2 = extractfield(pb_f1, 'iteration') <= x_max;
                pb_f2 = pb_f1(i_2);
                i_3 = extractfield(pb_f2, 'cons_identity') >= y_min;   
                pb_f3 = pb_f2(i_3);

                % Find 'coi' values.
                align_length = extractfield(pb_f3, 'align_homo_align_length');
                values = extractfield(pb_f3, coi_name);
                norm = values./align_length;

                % Generate container arrays for all experiments.
                COI = [COI; norm'];
                
            end
        end
        
        % Determine number of total events in this data set.
        SE(e) = length(COI);

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
    
    end

    % Build grouping variable array.
    GE = []; 

    for j = 1:length(list)

        % Initialize containers.
        ge = [];

        for k = 1:SE(j)
            ge{k} = j;
        end

        % Update concatinated grouping variable array.
        GE = [GE; ge'];

    end

    % Create dwell time (s) boxplot of all filtered events (exp).
    boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b');
    %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1); 
    grid;
    title(['filtered ', ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name, '_filter'];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end    
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Kinetic alignment filter end');
fprintf('\n');

% End timer.
toc