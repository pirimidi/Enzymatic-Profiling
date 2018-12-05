function [barcodes, scores] = barcode_check_homo(barcodes, ran_bar, threshold, scores)

% Set up regular expressions for "lumping" "stuttery" regions.
expr = {'A+', 'C+', 'T+', 'G+'};
repl = {'A'; 'C'; 'T'; 'G'};

% Universal primer sections flanking the random barcode section (middle).
l = 'GGCTAAAAT';
r = 'TCCCCACTCT';

% Initialize parameters.
score = 0; alignment = 0; identity = 100;

% Initialize a condition counter to be a non-empty vector.
filter = [0];
    
% Determine the barcode comparison pairs.
v = 1:length(barcodes);
C = nchoosek(v, 2); % all comparisons in current set
I = find(C(:,2) == length(barcodes)); % find comparisons revevant to current level
    
while (isempty(filter) == 0)    
  
    % Update new barcode.
    barcode = [l randseq(ran_bar) r];
        
    % Perform all barcode-to-barcode alignment comparison, which involves 
    % the newly generated barcode.
    for i = 1:length(I)

        % Current comparison pair.
        c = C(I(i), :);
        j = c(1); k = c(2);
        
        lumped_template_j = regexprep(barcodes{j}, expr, repl, 'ignorecase'); % barcode = regexprep(raw_template, expr, repl, 'ignorecase'); % barcode
        lumped_template_k = regexprep(barcode, expr, repl, 'ignorecase'); % barcode = regexprep(raw_template, expr, repl, 'ignorecase'); % barcode

        [score, alignment] = swalign(lumped_template_j, lumped_template_k);   
        total = length(alignment(2, :));
        match = length(find(alignment(2, :) == '|'));
        identity = match/total*100;  
        
        % Store alignment identity score value for the comparison.
        scores(j,k) = identity;
        
    end

    % Check if all barcode alignment comparisons meet the threshold crietria.
    filter = scores(scores > threshold);

end

% Update barcode set with new element.
barcodes{k} = barcode;
