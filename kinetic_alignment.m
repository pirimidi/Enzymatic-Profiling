%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: March 6, 2017
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: This program receives a set of MAT-files, one for each ('pol-bar') 
% experiment, containing all parameters parsed from the 'cell_annotations' 
% EXCEL file. Then, for each 'pol-bar', it iterates through the experimental 
% folders ('16XXXX_HMS_XX_xxxxxxx_XXXXXXXXX') and generates normalized
% whisker plots (4-4 each) of alignment accuracy for:
%
%  - heteropolymer deletions/identicals/isertions/mismatches
%  - homopolymer deletions/identicals/isertions/mismatches
%
% DEFINITIONS:      
%
% (1) hetero alignment accuracy (%) = align_num_ident / align_align_length
% (2) homo alignment accuracy (%) = align_homo_num_ident / align_homo_align_length
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_alignment

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic alignment start');
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
data_dir = strcat(work_dir, '\data\mat_files\stats');
cd(data_dir);

disp(['--> IN DIRECTORY: ' data_dir]);

% Read in all 'experiment' folder names one-by-one.
list = dir('pol6*_stats*.mat');

% Figure counter.
counter = 1;

disp('--> NORMALIZED WHISKER PLOTTING SECTION');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              WHISKER PLOTTING                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate 'experimental set' directory to hold the whisker plots.
cdir = 'kinetic_alignment';

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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                    
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
    savefig([fn, '.fig']);
    print('-dbmp', [fn, '.bmp']);
    
    % Update counter.
    counter = counter + 1;

    disp(['--> PROCESSED DATA SET: ' fn]);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NNormalized whisker plot for heteropolymer identicals.
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_homo_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_homo_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_homo_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
        for exp = 1:length(cells_list(bar).experiment)

            % Only look at structures with elements (can be empty).
            if length(cells_list(bar).experiment(exp).s) > 1
                
                % Iterate through all experiments per (pol-)bar and plot desired field. 
                pb = cells_list(bar).experiment(exp).s;

                % Find 'coi' values.
                align_length = extractfield(pb, 'align_homo_align_length');
                values = extractfield(pb, coi_name);
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
    title([ordinate_name, ' bar-kinetics']);
    xlabel('pol id');
    ylabel(ordinate_name);
    ylim([0 1]);

    % Save plot.
    fn = [num2str(counter), '_', barcodes{bar}, '_', coi_name];
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
disp('--> Kinetic alignment end');
fprintf('\n');

% End timer.
toc