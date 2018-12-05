function [temp] = random_barcode(ran_bar)

temp = randseq(ran_bar); % random 32-mer 

% Generate sequence without homopolymer runs.
for n = 1:ran_bar-1        
    if temp(n) == 'A' 
        alphabet = {'C', 'T', 'G'};            
    elseif temp(n) == 'C' 
        alphabet = {'A', 'T', 'G'};            
    elseif temp(n) == 'T' 
        alphabet = {'A', 'C', 'G'};            
    else
        alphabet = {'A', 'C', 'T'};
    end

    temp(n+1) = alphabet{randi([1 3], 1)}; % update base at query position
end

            