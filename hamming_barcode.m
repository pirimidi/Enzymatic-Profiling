function [temp] = hamming_barcode(barcode_j, barcode_k, ran_bar)

disp('--> In Hamming barcode!');

% Set up regular expressions for "lumping" "stuttery" regions.
expr = {'A+', 'C+', 'T+', 'G+'};
repl = {'A'; 'C'; 'T'; 'G'};

% Lump barcodes and modify match positions to meet treshold. 
lumped_barcode_j = regexprep(barcode_j, expr, repl, 'ignorecase');
lumped_barcode_k = regexprep(barcode_k, expr, repl, 'ignorecase');

[score, alignment, start] = nwalign(lumped_barcode_j, lumped_barcode_k);  
total = length(alignment(2, :));
match = length(find(alignment(2, :) == '|'));
identity = match/total*100; 

start(1)
alignment
length(alignment)

% Initialize match counters and index array.
mat_count = 0; del_count = 0;
mat_index = []; del_index = [];

% Compute number of matches.
for k = 1:length(alignment)
   
    % Count number of matches.
    if strcmp(alignment(1+3*(k-1)), alignment(3+3*(k-1)))
        
        mat_count = mat_count + 1;
        mat_index(mat_count) = k;
        
    % Count number of deletions.    
    elseif strcmp(alignment(3+3*(k-1)), '-')

        del_count = del_count + 1;
        del_index = [del_index k];
        
    end
    
end

aligned_barcode_k = alignment(3+3*(0:length(alignment)-1));

% Replace all matches with a random mismatch.
for n = 1:length(mat_index)     
    if 6 < mat_index(n) && mat_index(n) < length(alignment)-7

        base = alignment(3+3*(mat_index(n)-1));

        if base == 'A' 
            alphabet = {'C', 'T', 'G'};            
        elseif base == 'C' 
            alphabet = {'A', 'T', 'G'};            
        elseif base == 'T' 
            alphabet = {'A', 'C', 'G'};            
        else
            alphabet = {'A', 'C', 'T'};
        end
        
        aligned_barcode_k(mat_index(n)) = alphabet{randi([1 3], 1)}; % update base at query position

    end
end

% Replace all deletions with a random mismatch.
for n = 1:length(del_index)     
    if 6 < del_index(n) && del_index(n) < length(alignment)-7
        
        base = alignment(1+3*(del_index(n)-1));

        if base == 'A' 
            alphabet = {'C', 'T', 'G'};            
        elseif base == 'C' 
            alphabet = {'A', 'T', 'G'};            
        elseif base == 'T' 
            alphabet = {'A', 'C', 'G'};            
        else
            alphabet = {'A', 'C', 'T'};
        end
        
        aligned_barcode_k(del_index(n)) = alphabet{randi([1 3], 1)}; % update base at query position

    end
end

% Remove gaps.
unlumped_barcode_k = regexprep(aligned_barcode_k, '-', '');

% Return only the lumped Hamming edited barcode (middle).
lumped_temp = unlumped_barcode_k((6:length(lumped_barcode_k)-7));

% Add back homopolymer region to lumped barcode end. (Will be removed with
% lumping eventually, so does not make any difference.)
if(length(lumped_temp) < ran_bar)
    
    last_base = lumped_temp(end);
    iter = ran_bar - length(lumped_temp);
    for i = 1:iter
        lumped_temp = [lumped_temp last_base];
    end
    
end

% Return unlumped 32 bp barcode region.
temp = lumped_temp;
