%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: June 13, 2017.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a standard barcode design: 
% 5'-GGCTAAAAT-[32-bp barcode]-TCCCCACTCT-3', this is program generates 96 
% unique templates by seeding from an initial template and recursively 
% compares them by alignment for uniqueness. A new random template is 
% generated if a particular alignment has greater than 50% base identity, 
% otherwise it is accepted as a distinguishable template-barcode. This 
% iterative barcode generation and alignment procedure continues until all 
% template-barcodes are unique.
%
% Input arguments: 'num_bar' = number of desired barcodes
%                  'ran_bar' = number of random bases in the template 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function recursive_naive_barcodes_nw(num_bar, ran_bar, threshold)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Recursive naive barcodes NW start');
fprintf('\n');

% Set default number formatting.
format short;

% Turn off warnings during run.
warning('off', 'all');

% Define current working directory.
work_dir = pwd;

%-------------------------------------------------------------------------%
%                         RANDOM BARCODE GENERATION                       %
%-------------------------------------------------------------------------%

% Define barcode cell array to hold randomly generated barcodes.
barcodes = {};

% Generate a barcode alignment identity score matrix for each comparison.
scores = zeros(num_bar);

% Universal primer sections flanking the random barcode section (middle).
l = 'GGCTAAAAT';
r = 'TCCCCACTCT';

disp('--> RANDOM BARCODE GENERATION SECTION');
 
% Generate initial random (seed) barcode according to template-barcode design.
temp = random_barcode(ran_bar);
barcodes{1} = [l temp r];

%-------------------------------------------------------------------------%
%                        SEQUENCE ALIGNMENT OUTPUT                        %
%-------------------------------------------------------------------------%

while(length(barcodes) ~= num_bar)
    
    disp(['--> Calculating barcode: ', num2str(length(barcodes))]);
    temp = random_barcode(ran_bar);
    barcodes{length(barcodes)+1} = [l temp r];
    [barcodes, scores] = barcode_check_nw(barcodes, ran_bar, threshold, scores);
    
end

scores

% Unique barcode data structure.
bc_data = {};

% Generate FASTA file for unique barcodes.
for m = 1:length(barcodes)
    
    bc_data(m).Header = num2str(m);
    bc_data(m).Sequence = barcodes{m};
    
end

% Generate FASTA file for barcode comparison history..
fastawrite('unique_barcodes_nw.fasta', bc_data);

disp('--> Recursive naive barcodes NW end');
fprintf('\n');

% End timer.
toc
