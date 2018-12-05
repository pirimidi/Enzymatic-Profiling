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
% whisker plots featuring alignment accuracy (%) versus number of full 
% barcode iteration ("align_copies") for each 'pol-bar' case displaying 
% them on the same plot for comparison. NOTE: Use MATLAB R2017a to run this
% code, older version might trigger errors. 
%
% DEFINITIONS:      
%
% (1) hetero alignment accuracy (%) = align_num_ident / align_align_length
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function kinetic_cutoff_all

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Kinetic cutoff polymerase fix start');
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
data_dir = strcat(work_dir, '/data/mat_files/stats');
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
cdir = 'kinetic_cutoff_all';

% Define direcory to hold figures.
if ~exist(cdir, 'dir')
  mkdir(cdir);
end

% Iterate through all 'pol' folders.
 for e = 1:length(list) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalized whisker plot for heteropolymer identicals.
    coi_name = 'align_num_ident';
    ordinate_name = regexprep(coi_name, '_', '-');
    
    % Define array container for fixed 'bar' and grouping variable array.
    BAR = []; GE = [];

    % Iterate through all three barcodes.
    for bar = 1:length(barcodes)  

        % Load nanopore data structure per experimental set (pol-bar).
        load(list(e).name);
        disp(['--> MAT-FILE LOADED: ' list(e).name]);
        ss = strsplit(list(e).name, '.');
        
        % Define array container for 'coi' and grouping variable array.
        COI = []; ge = [];

        % Iterate through all 'experiment' folders for (pol-)bar.
        for exp = 1:length(cells_list(bar).experiment)

            % Iterate through all experiments per (pol-)bar and plot desired field. 
            pb = cells_list(bar).experiment(exp).s;

            % Find align copy values for each 'coi'.
            align_copies = extractfield(pb, 'align_copies');
            
            % Update concatinated grouping variable array.
            ge = [ge; align_copies'];
            
            % Find 'coi' values.
            align_length = extractfield(pb, 'align_align_length');
            values = extractfield(pb, coi_name);
            norm = values./align_length;

            % Generate container arrays for all experiments.
            COI = [COI; norm'];
            
        end

        % Generate container array for 'pol-bar' combination (fix bar).
        BAR = [BAR; COI];
        
        % Update concatinated grouping variable array.
        GE = [GE; ge];
        
        % Generate accumulative plot. 
        
        
        % Generate figure for each pol-bar combination.
        figure(counter);
    
        % Create dwell time (s) boxplot of all filtered events (exp).
        boxplot(BAR, GE, 'PlotStyle', 'compact', 'Symbol', 'b', 'Whisker', 0.2413);
        %boxplot(BAR, GE, 'PlotStyle', 'compact', 'OutlierSize', 1, 'Whisker', 0.2413); 
        grid;
        title([ordinate_name, ' iteration-cutoff-all HETERO']);
        xlabel('number of iterations');
        ylabel(ordinate_name);
        ylim([0 1]);

        % Save plot.
        pp = strsplit(ss{1}, '_');
        fn = [pp{1}, '_', barcodes{bar}, '_', coi_name, '_HETERO'];
        savefig([fn, '.fig']);
        print('-dbmp', [fn, '.bmp']);

        % Update counter.
        counter = counter + 1;

        % Summary statistics organized by group.
        [mean, sem, numel, gname, std, var, min, max, range] = ...
         grpstats(BAR, GE, {'mean', 'sem', 'numel', 'gname', 'std', 'var', ...
                            'min', 'max', 'range'})

        disp(['--> PROCESSED DATA SET: ' fn]); 
    
    end 
end    
    
% Move all figures to 'plots' directory.
movefile('*.fig', cdir);
movefile('*.bmp', cdir);

% Close all open figures.
close all;        

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Kinetic cutoff all end');
fprintf('\n');

% End timer.
toc