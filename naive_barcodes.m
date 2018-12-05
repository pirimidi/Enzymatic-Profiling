%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: May 1, 2017.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a standard barcode design: 
% 5'-GGCTAAAAT-[32-bp barcode]-TCCCCACTCT-3', this is program generates 96
% random templates then iterates through all of them and compares them by
% alignment for uniqueness. A new random template is generated if a particular 
% alignment has greater than 60% base identity, otherwise it is accepted as
% a distinguishable template-barcode. This iterative barcode generation and 
% alignment procedure continues until all template-barcodes are unique.
%
% Input arguments: 'num_bar' = number of desired barcodes
%                  'ran_bar' = number of random bases in the template 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function naive_barcodes(num_bar, ran_bar, threshold)

%-------------------------------------------------------------------------%
%                                 STARTUP                                 %
%-------------------------------------------------------------------------%

% Start timer.
tic 

fprintf('\n');
disp('--> Naive barcodes start');
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

% Universal primer sections flanking the random barcode section (middle).
l = 'GGCTAAAAT';
r = 'TCCCCACTCT';

disp('--> RANDOM BARCODE GENERATION SECTION');
 
% Generate 96 random barcode according to template-barcode design.
for i = 1:num_bar
    barcodes{i} = [l randseq(ran_bar) r];
end

% Generate a barcode alignment identity score matrix for each comparison.
scores = zeros(num_bar);

% Initialize a condition counter to be a non-empty vector.
filter = [0];

%-------------------------------------------------------------------------%
%                        SEQUENCE ALIGNMENT OUTPUT                        %
%-------------------------------------------------------------------------%

% Alignment iterator.
ai = 1;

% Alignment comparison data structure.
ac_data = {};

% Initialize parameters,
score = 0; alignment = 0; identity = 100;

while (isempty(filter) == 0)    
    for j = 1:num_bar    
        
        % Compare current (j) template to all other templates in the queue
        % below (n<j) by alignment.
        for k = 1:num_bar 
            
            % Do alignment comparison if the current barcode index is less than
            % the other barcode in the queue.
            if j < k
                while (identity > threshold)
                    
                    barcodes{k} = [l randseq(ran_bar) r];
                    [score, alignment] = swalign(barcodes{j}, barcodes{k});   
                    disp(['--> Raw score ' num2str(j), '-', num2str(k), ': ' num2str(score)]);
                    total = length(alignment(2, :));
                    match = length(find(alignment(2, :) == '|'));
                    identity = match/total*100;
                    disp(['--> Identity: ' num2str(identity) '%']);
                    alignment  

                    % Store alignment identity score value for the comparison.
                    scores(j,k) = identity;                    
                end
                
                % Write alignment details into FASTA format.
                ac_data(ai).Header = [num2str(ai), ': ' num2str(j), '-', num2str(k), ', identity: ', num2str(identity), '%, template: ', barcodes{j}];
                ac_data(ai).Sequence = barcodes{k};
                ac_data(ai).Template = barcodes{j};
                ac_data(ai).Score = score;
                ac_data(ai).Identity = identity;
                ac_data(ai).ID = [num2str(j), '-', num2str(k)]; 
                
                % Increment counter.
                ai = ai + 1;

                % Reset identity parameter,
                identity = 100;            
            end  
        end
    end

    % Check if all barcode alignment comparisons meet the threshold crietria.
    filter = scores(scores > threshold);
    
end

scores

% Generate FASTA file for barcode comparison history..
fastawrite('unique_comparisons.txt', ac_data);

% Unique barcode data structure.
bc_data = {};

% Generate FASTA file for unique barcodes.
for m = 1:length(barcodes)
    
    bc_data(m).Header = num2str(m);
    bc_data(m).Sequence = barcodes{m};
    
end

% Generate FASTA file for barcode comparison history..
fastawrite('unique_barcodes.fasta', bc_data);

disp('--> Naive barcodes end');
fprintf('\n');

% End timer.
toc
