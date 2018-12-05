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

function delete_files(folder_name)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Delete files start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working and data directory.
work_dir = pwd;

% Define barcode list.
barcodes = {'comp3', 'fv2', 'rep3'};
%barcodes = {'fv2'};

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
        
        % Change directory to experimental sub-folder.
        cd([data_dir, '\', list(exp).name, '\plots']);
        
        % REMOVE after using it once!
        delete('*');
        
        % Move up a directory.
        cd('../');
    end   
end

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Delete files end');
fprintf('\n');

% End timer.
toc