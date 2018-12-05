%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 31, 2015.
%
% For: Climathon scoring algorithm with Claudine's team.
%
% Purpose: This program receives a parcel ID (or GPS location) and finds
% all building characteristic values (d_i,j) in the 'data matrix' 
% associated with scoring categories (sc) tabulated in the 'scoring matrix' 
% (prepared by Claudine). Then, using these parameters, it finds the 
% associated weight (w_i,j)and relative score (r_i,j) for each in the 
% 'scoreing matrix' and calculates the final - weighted - score according 
% to formula:
%
%              ws = sum(w_i,j * r_i,j), j = 1 to max(sc)
%
% ..., where 'i' is the parcel ID (row identifier), while 'j' is the list 
% of categories used for scoring.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function scoring_algorithm_m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            SCORING STARTUP                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start timer.
tic

fprintf('\n');
disp('--> Scoring algorithm start');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

% Define direcory to hold figures.
if ~exist('plots', 'dir')
  mkdir('plots');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         SCORING MATRIX PARSING                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read data scoring categories from 'scoring matrix'.
[NM, TX, RW] = xlsread('scoring_matrix.xlsx');

RW

% Create an empty cell array for category keyset/valueset in hash table.
category = {}; value_set = {};

% Create a data structure to store weight/relative score values.
values = struct;

% Initialize category counter.
c = 1;

% Parse 'scoring matrix' and create a data structure with its contents.
for i = 1:length(RW)
    
    % If the weight is a 'NaN', it is a category.
    if isnan(RW{i,2}) == 1
        disp(['--> Processing category: ', RW{i,1}]);
        
        % Initialize value counter.
        v = 1;
        
        % If we passed the first category, assign values structure to value
        % set in hash table.
        if c-1 ~= 0
            
            values.value = value;
            values.weight = weight;
            values.r_score = r_score;
            
            value_set{c-1} = values;
        end
        
        % Initialize value/weight/relative score arrays.
        value = {}; weight = []; r_score = [];
        
        % Add to category cell array.
        category{c} = RW{i,1};
        
        % Increment category counter.
        c = c+1;
        
    else
        
        % Generate arrays to hold value/weight/relative score pairs for the
        % particular category.
        value{v} = RW{i,1};
        weight = [weight RW{i,2}];
        r_score = [r_score RW{i,3}];
        
        % Increment value counter.
        v = v+1;
        
    end
end
   
% If we passed the last category, assign the last values structure to value
% set in hash table.
values.value = value;
values.weight = weight;
values.r_score = r_score;

value_set{c-1} = values;

% Build a hash table. 
map = containers.Map(category, value_set);

% Line separator.
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         SCORING MATRIX PARSING                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in data scoring categories from 'data matrix'.
[ND, TD, RD] = xlsread('data_matrix.xlsx');

[r,c] = size(RD);

% Generate category list from 'data matrix'.
for j = 1:c
    CA{j} = RD{1,j};
end

% Generate parcel id list from first column of 'data matrix'.
for k = 1:r
    PA{k} = RD{k,1};
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            SCORE CALCULATION                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Create a data structure to store partial/final score values.
scores = struct;

% Iterate through all parcel ids and calculate partial/final score. 
for q = 2:length(PA)
    
    % Current parcel id. 
    parcel_id = PA{q};
    
    % Find row of parcel id.
    id_p = find(strcmp(parcel_id, PA'));

    % Match up 'data matrix' query with 'scoreing matrix' values.
    categories = keys(map);

    % Initialize weight and relative score arrays for partial scores.
    weight_f = []; r_score_f =[];

    for l = 1:length(categories)

        % Look up value and associated weight/relative score by category.
        v = map(categories{l});

        % Initialize cell array for value list of categories.
        sub_value = {};

        % Find column of scoring category.
        id_c = find(strcmp(categories{l}, CA'));
        
        % Find value in 'data matrix' corresponding to pardel id/category pair.
        data_value = RD{id_p, id_c};

        % Data matrix entry associated with parcel id/category pair is zero.
        if data_value == 0

            % Assign zeros for weight/score in this case.
            weight_c = 0;
            r_score_c = 0;

        % Data matrix entry associated with parcel id/category pair is 'NaN'.
        elseif isnan(data_value) == 1

            % Assign zeros for weight/score in this case.
            weight_c = 0;
            r_score_c = 0;

        % If there is a valid data matrix entry associated with parcel
        % id/category pair.
        else 

            % Regenerate category value list, depending on type.
            if strcmp(categories{l}, 'R_EXT_FIN') || strcmp(categories{l}, 'R_HEAT_TYP') ...
                                                  || strcmp(categories{l}, 'S_EXT_FIN') 
                for m = 1:length(v.value)   
                    sub_value{m} = v.value{m}(1);
                end

            else

                % If not in special category, proceed with simple match case.
                sub_value = v.value;

            end

            % Find matching value in 'scoring matrix' with associated
            % weight/relative score.
            if strcmp(categories{l}, 'YR_BUILT') || strcmp(categories{l}, 'YR_REMOD')

                for n = 1:length(sub_value)

                    % Get range start/end values for each sub-value.
                    range = str2num(sub_value{n});

                    % Check which range is the data value in.
                    if range(1) <= data_value && data_value <= range(2)
                        id_s = n;
                        break
                    end

                end

            % Simple string match case.
            else
                id_s = find(strcmp(data_value, sub_value));
            end

            weight_c = v.weight(id_s);
            r_score_c = v.r_score(id_s); 

        end
        
        % And the category score is...
        scores.category_score(q, l) = weight_c * r_score_c;

        weight_f = [weight_f weight_c];
        r_score_f = [r_score_f r_score_c];

    end

    % And the final score is...
    disp(['--> Processing parcel id: ' parcel_id]);
    scores.final_score(q) = weight_f * r_score_f';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCORE INPUT                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Match up 'data matrix' query with 'scoreing matrix' values.
categories = keys(map);

for b = 1:length(categories)

    % Look up value and associated weight/relative score by category.
    v = map(categories{b});

    % Determine maximum weigth/relative scores for each category.
    max_weigth(b) = max(v.weight);
    max_r_score(b) = max(v.r_score);

end

% Theoretical maximum score.
max_score = max_weigth * max_r_score';

% Normalize score to 1.
norm_final_score = scores.final_score/max_score;

% Scale score to [-5,+5] according to (y2-y1)/(x2-x1)*x + b - equation of line.
resili_score = 10*norm_final_score-5;

% Generate new columns in data matrix with category and final scores.
xlswrite('data_matrix.xlsx', resili_score(2:end)', 'TaxFloodParcels', 'AW2');
 
% Generate new columns in data matrix with category and final scores.
xlswrite('data_matrix.xlsx', scores.final_score(2:end)', 'TaxFloodParcels', 'AX2');
 
flood_score = scores.category_score(2:end, 1);
xlswrite('data_matrix.xlsx', flood_score, 'TaxFloodParcels', 'AY2');
 
r_ext_fin_score = scores.category_score(2:end, 2);
xlswrite('data_matrix.xlsx', r_ext_fin_score, 'TaxFloodParcels', 'AZ2');
 
r_heat_type_score = scores.category_score(2:end, 3);
xlswrite('data_matrix.xlsx', r_heat_type_score, 'TaxFloodParcels', 'BA2');
 
yr_built_score = scores.category_score(2:end, 4);
xlswrite('data_matrix.xlsx', yr_built_score, 'TaxFloodParcels', 'BB2');
 
yr_remod_score = scores.category_score(2:end, 5);
xlswrite('data_matrix.xlsx', yr_remod_score, 'TaxFloodParcels', 'BC2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              SCORE PLOTTING                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot rate constants as a function of starting mass.
x = 1:length(resili_score(2:end)); 

% ResiliScore plot.
figure(1);
plot(x, resili_score(2:end)', 'rx', 'MarkerFaceColor', 'r', 'MarkerSize', 2);
grid;
title('ResiliScore');
xlabel('Parcel IDs');
ylabel('Score');
axis([-500 6800 -7 6]);
legend('ResiliScore [-5,+5]');

% Save plot.
savefig('resili_score.fig');
print('-dbmp', 'resili_score.bmp');

% Unscaled composite plot.
figure(2);
plot(x, scores.final_score(2:end)', 'bx', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
hold on;
plot(x, flood_score, 'gx', 'MarkerFaceColor', 'g', 'MarkerSize', 2);
hold on;
plot(x, r_ext_fin_score, 'mx', 'MarkerFaceColor', 'm', 'MarkerSize', 2);
hold on;
plot(x, r_heat_type_score, 'cx', 'MarkerFaceColor', 'c', 'MarkerSize', 2);
hold on;
plot(x, yr_built_score, 'yx', 'MarkerFaceColor', 'y', 'MarkerSize', 2);
hold on;
plot(x, yr_remod_score, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 2);

% Decorate plot.
grid;
title('Unscaled Final/Category Score(s)');
xlabel('Parcel IDs');
ylabel('Score');
axis([-500 6800 -0.5 5.5]);

legend('Final Score [0,5]', 'Flood Score', 'Exterior Finish Score', ...
       'Heating Type Score', 'Year Built Score', ...
       'Year Remodeled Score', 'Location', 'northeastoutside');

% Save plot.
savefig('unscaled_score.fig');
print('-dbmp', 'unscaled_score.bmp')

% Change plot background to light grey.
whitebg([0.9 0.9 0.9]);

% Move all figures to 'plots' directory.
movefile('*.fig', 'plots');
movefile('*.bmp', 'plots');

fprintf('\n');
disp('--> Scoring algorithm end');
fprintf('\n');

% End timer.
toc
