function [Result,AdaptFuncValue] = ASMRDE(PopSize,SearchDimension,SScope,FuncNo,MaxGen)

% Adaptive Social Mobility-Restructuring Differential Evolution 
% Compiler: Yiwen Zhuo

% Result: best solution
% AdaptFuncValue: Fitness change procee
    
% Define random seed
rand('state',sum(100*clock));
% For CEC test, PopSize_Input equals problem dimension
sample_size = SearchDimension;
% Define maximum function evaluations
NFEmax = SearchDimension*10000;
% Define initial population
POPS = rand(PopSize,SearchDimension+1);
% Random population initialization
for SearchDimension = 1:SearchDimension
    if SScope(SearchDimension,1)==-inf && SScope(SearchDimension,2)~=inf
        POPS(:,SearchDimension) = unifrnd(SScope(SearchDimension,2)-10^10,SScope(SearchDimension,2),PopSize,1);
    end

    if SScope(SearchDimension,1)~=-inf && SScope(SearchDimension,2)==inf
        POPS(:,SearchDimension) = unifrnd(SScope(SearchDimension,1),SScope(SearchDimension,1)+10^10,PopSize,1);
    end

    if SScope(SearchDimension,1)==-inf && SScope(SearchDimension,2)==inf
        POPS(:,SearchDimension) = unifrnd(-10^10,10^10,PopSize,1);
    end

    if SScope(SearchDimension,1)~=-inf && SScope(SearchDimension,2)~=inf
        POPS(:,SearchDimension) = POPS(:,SearchDimension)*(SScope(SearchDimension,2)-SScope(SearchDimension,1))+SScope(SearchDimension,1);
    end
end
% For CEC test problems, FuncNo is the function number
POPS(1:PopSize,SearchDimension+1) = cec17_func(POPS(1:PopSize,1:SearchDimension)',FuncNo)';
% Find individual with the best fitness value
[~,row] = min(POPS(1:PopSize,SearchDimension+1));
% Define BEST as the best individual
BEST = POPS(row,1:SearchDimension+1);
% Record fitness values of evaluated individuals for later iteration curve plotting
nn = 1;
SearchProcess(1,nn) = POPS(1,SearchDimension+1);
nn = nn + 1;

for ii = 2:PopSize
    if SearchProcess(1,nn-1) > POPS(ii,SearchDimension+1)
        SearchProcess(1,nn) = POPS(ii,SearchDimension+1);
        nn = nn + 1;
    else
        SearchProcess(1,nn) = SearchProcess(1,nn-1);
        nn = nn + 1;
    end
end
% Initial evaluation count
nfes = PopSize;

% Define initial iteration count
k = 1;

% Range of parameter F variation
delta = 0.025;

% Initialize archive B
B = [];
POPM = [POPS(1:end,:);B];

% Set maximum number of identical individuals
max_duplicates = floor(PopSize / 4);

% Set cycle period
Loop = 15;

% Diversity threshold
div_threshold = 0.74;

while nfes < NFEmax

    % Update iteration count
    k = k + 1;
    
    % Update evaluation count
    nfes = nfes + PopSize;

    % Update temperature
    T = 17.8 * (1 - k / MaxGen) + 0.2; % Linear decrease
    
    % Social Restructuring once every Loop iterations
    if mod(k,Loop) ==0 && k ~= MaxGen 
        % Selection operation
        uniquePOPM = unique(POPM(:,:), 'rows', 'stable');
        selection_prob_unique = calculateSelectionProbability(uniquePOPM, size(uniquePOPM, 1), SearchDimension, T);
        [new_POPS, new_POPS_indices, B] = selectAndUpdatePopulation_uniquePOPM(POPM, B, PopSize, SearchDimension, selection_prob_unique);
    end
        % Establish population union P¡ÈB
        POPM = [POPS;B];
        Diversity_record_B(k-1) = calculateDiversity(B);
        Diversity_record_POPM(k-1) = calculateDiversity(POPM);

        % Evolution following population best
        [uniquePOPS, ~, ~] = unique(POPS(:,:), 'rows', 'stable');

        uniquePOPSsolution = size(uniquePOPS,1);

        % Select Xpbest base vector, first calculate parameter p: elite ratio, with linear decreasing change
        [~,indexbest] = sort(uniquePOPS(1:end,1+SearchDimension),'ascend');
        p = 0.4 - 0.2*(nfes/NFEmax)^1;
        pNP = round(max(1,uniquePOPSsolution*p));

        % Generate random numbers to select Xpbest
        randindex = ceil(rand(PopSize,1)*pNP);
        randindex = indexbest(randindex);
        Xpbest = uniquePOPS(randindex,1:SearchDimension);
      
        % Calculate F value
        MeanF = 0.7 - (0.7 - 0.025) * (nfes-1)/(NFEmax-1);
        POPSdiversity = calculateDiversity(POPS(1:PopSize,1:SearchDimension));
        if POPSdiversity >=div_threshold*PopSize
            MeanF = max(MeanF,0.5);
        end

        F = normrnd(MeanF, delta/2, PopSize, 1);
        % Ensure F values are within reasonable range (e.g., between 0 and 1)
        F = max(0, min(1, F));

        if POPSdiversity >= div_threshold*PopSize && k >= 0.5*MaxGen
            meanCR = 0.92;
            CR = normrnd(meanCR, 0.1, PopSize, 1);
        else
            meanCR = 0.9;
            CR = normrnd(meanCR, 0.1, PopSize, 1);
        end
        
        CR = max(0, min(1, CR));
        
        % Record diversity
        Diversity_record(k-1) = POPSdiversity;
        % Record MeanF
        MeanF_record(k-1) = MeanF;
        
    if mod(k,Loop) ==0 && k ~= MaxGen 
        rr = selectRandomIndices3BasedOnNewPOPS(PopSize, new_POPS_indices, POPM);   
        V = generateMutationVector(new_POPS, POPM, Xpbest, F, SScope, PopSize, SearchDimension, rr);
        POPS = new_POPS;
    else        
        rr = selectRandomIndices3(PopSize,POPM);  
        V = generateMutationVector(POPS, POPM, Xpbest, F, SScope, PopSize, SearchDimension, rr);    
    end

    % Define trial vector U
    U = POPS(1:PopSize,1:SearchDimension+1);

    % Randomly select one dimension that must be crossed
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);

    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);

    % I is a logical matrix
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);

    % Evaluate trial vectors
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',FuncNo);

    tmp = (U(1:PopSize,SearchDimension+1) < POPS(1:PopSize,SearchDimension+1));        % Generate tmp as a marker for offspring replacing parents

    temp = repmat(tmp,1,SearchDimension+1);                                                % D-dimensional copy version of tmp

    % Update population
    POPS(1:PopSize,1:SearchDimension+1) = temp.*U(1:PopSize,1:SearchDimension+1) + (1-temp).*POPS(1:PopSize,1:SearchDimension+1);

    % Record function evaluation values
    for i = 1:PopSize
        SearchProcess(1,nn) = ...
            (SearchProcess(1,nn-1) > POPS(i,SearchDimension+1)).*POPS(i,SearchDimension+1) +...
            (SearchProcess(1,nn-1) <= POPS(i,SearchDimension+1)).*SearchProcess(1,nn-1);
        nn = nn + 1;
    end
    
    % Update merged population
    POPM = [POPS(1:end,:);B];

    % Select and record the best individual from the merged population
    [~,row] = min(POPM(1:end,SearchDimension+1));
    BEST = POPM(row,1:SearchDimension+1);

end

% Define generation interval for iteration curve plotting
kk = 1:10000;

% Output result AdaptFuncValue for curve plotting
AdaptFuncValue = SearchProcess(sample_size.*kk);

% Output result Result for calculating Mean, SD, Best, Worst values, used for statistical analysis
Result = BEST;

% Function section

% Calculate population diversity
function diversity = calculateDiversity(population)
    [uniquePopulation, ~, ~] = unique(population, 'rows');
    diversity = size(uniquePopulation, 1);
end

% Calculate selection probability
function selection_prob = calculateSelectionProbability(POPS, PopSize, SearchDimension, T)
    % Calculate fitness (directly use original fitness values, smaller is better)
    fitness = -POPS(1:PopSize, SearchDimension+1);

    % Normalize fitness (Z-score normalization), but handle cases where standard deviation is zero
    fitness_mean = mean(fitness);
    fitness_std = std(fitness);

    if fitness_std == 0
        % If standard deviation is zero, all fitness values are the same
        % In this case, we can give each individual equal selection probability
        selection_prob = ones(PopSize, 1) / PopSize;
    else
        % Normal normalization and subsequent calculations
        fitness_normalized = (fitness - fitness_mean) / fitness_std;
        
        % Use softmax to calculate selection probability
        exp_fitness = exp(fitness_normalized / T);
        selection_prob = exp_fitness / sum(exp_fitness);
        
        % Add small positive ¦Å to ensure no zero probability
        epsilon = 1e-10;
        selection_prob = selection_prob + epsilon;
        selection_prob = selection_prob / sum(selection_prob);
    end
end

function B = updateArchive(B, POPS, unselected_indices, SearchDimension, PopSize)
    % maxArchiveSize = ceil(1.5 * PopSize);
    maxArchiveSize = PopSize;
    % maxArchiveSize = 0;
    
    for i = 1:length(unselected_indices)
        individual = POPS(unselected_indices(i), :);
        
        % If B is empty, add individual directly
        if isempty(B)
            B = individual;
        else
            % Check if the same individual already exists in B
            is_duplicate = any(all(B(:, 1:SearchDimension) == individual(1:SearchDimension), 2));
            
            % If not a duplicate, add to B
            if ~is_duplicate
                B = [B; individual];
            end
        end
    end

    % If B size exceeds maxArchiveSize, keep the best maxArchiveSize individuals
    if size(B, 1) > maxArchiveSize
        [~, sorted_indices] = sort(B(:, end));
        B = B(sorted_indices(1:maxArchiveSize), :);
    end
end

% Select two random numbers, requiring them to be different from each other and from the current individual
function rr = selectRandomIndices3(PopSize,POPM)
    % Pre-calculate unique identifier for each individual
    [~, ~, unique_ids] = unique(POPM(:,1:SearchDimension), 'rows');
    
    % Pre-calculate indices corresponding to each unique identifier
    unique_id_indices = arrayfun(@(x) find(unique_ids == x), 1:max(unique_ids), 'UniformOutput', false);
    
    nrandI = 2;
    rr = zeros(PopSize, nrandI);
    
    for i = 1:PopSize
        current_id = unique_ids(i);
        
        % Find all unique identifiers different from the current individual
        available_ids = setdiff(1:max(unique_ids), current_id);
        
        if length(available_ids) >= nrandI
            % Randomly select nrandI different identifiers
            selected_ids = available_ids(randperm(length(available_ids), nrandI));
            
            % For each selected identifier, randomly select a corresponding index
            selected = zeros(1, nrandI);
            for j = 1:nrandI
                possible_indices = unique_id_indices{selected_ids(j)};
                selected(j) = possible_indices(randi(length(possible_indices)));
            end
        else
            % If there aren't enough different individuals, use all available different individuals
            selected = cellfun(@(x) x(randi(length(x))), unique_id_indices(available_ids));
            
            % Fill remaining positions with random indices (may duplicate)
            while length(selected) < nrandI
                selected = [selected, randi(PopSize)];
            end
        end
        
        rr(i,:) = selected;
    end
end

% Select two random numbers, requiring them to be different from each other and from the current individual
function rr = selectRandomIndices3BasedOnNewPOPS(PopSize, new_POPS_indices, POPM)
    % Pre-calculate unique identifier for each individual
    [~, ~, unique_ids] = unique(POPM(:,1:SearchDimension), 'rows');
    
    % Pre-calculate indices corresponding to each unique identifier
    unique_id_indices = arrayfun(@(x) find(unique_ids == x), 1:max(unique_ids), 'UniformOutput', false);
    
    nrandI = 2;
    rr = zeros(PopSize, nrandI);
    
    for i = 1:PopSize
        % Use new_POPS_indices to map to original population index
        original_index = new_POPS_indices(i);
        current_id = unique_ids(original_index);
        
        % Find all unique identifiers different from the current individual
        available_ids = setdiff(1:max(unique_ids), current_id);
        
        if length(available_ids) >= nrandI
            % Randomly select nrandI different identifiers
            selected_ids = available_ids(randperm(length(available_ids), nrandI));
            
            % For each selected identifier, randomly select a corresponding index
            selected = zeros(1, nrandI);
            for j = 1:nrandI
                possible_indices = unique_id_indices{selected_ids(j)};
                selected(j) = possible_indices(randi(length(possible_indices)));
            end
        else
            % If there aren't enough different individuals, use all available different individuals
            selected = cellfun(@(x) x(randi(length(x))), unique_id_indices(available_ids));
            
            % Fill remaining positions with random indices (may duplicate)
            while length(selected) < nrandI
                selected = [selected, randi(size(POPM, 1))];
            end
        end
        
        rr(i,:) = selected;
    end
end

function V = generateMutationVector(new_POPS, POPM, Xpbest, F, SScope, PopSize, SearchDimension, rr)
    % Initialize mutation vector V
    V = zeros(PopSize,SearchDimension);

    % Generate V using formula DE/current-to-pbest/1 with external archive B
    % V(1:PopSize,1:SearchDimension) = new_POPS(1:PopSize,1:SearchDimension) + F * (Xpbest - new_POPS(1:PopSize,1:SearchDimension)) +...
    %     F*(POPM(rr(1:PopSize,1),1:SearchDimension) - POPM(rr(1:PopSize,2),1:SearchDimension));

    V(1:PopSize,1:SearchDimension) = new_POPS(1:PopSize,1:SearchDimension) + repmat(F(1:PopSize,1),[1,SearchDimension]).*(Xpbest - new_POPS(1:PopSize,1:SearchDimension)) +...
    repmat(F(1:PopSize,1),[1,SearchDimension]).*(POPM(rr(1:PopSize,1),1:SearchDimension) - POPM(rr(1:PopSize,2),1:SearchDimension));


    % Process dimensions in V that exceed boundaries, by taking the average of corresponding X and boundary range
    V(1:PopSize,1:SearchDimension) = ...
        ((V(1:PopSize,1:SearchDimension)>=repmat(SScope(:,1)',[PopSize,1]))&(V(1:PopSize,1:SearchDimension)<=repmat(SScope(:,2)',[PopSize,1]))).*(V(1:PopSize,1:SearchDimension))+...
        (V(1:PopSize,1:SearchDimension)<repmat(SScope(:,1)',[PopSize,1])).*((repmat(SScope(:,1)',[PopSize,1])+new_POPS(1:PopSize,1:SearchDimension))./2)+...
        (V(1:PopSize,1:SearchDimension)>repmat(SScope(:,2)',[PopSize,1])).*((repmat(SScope(:,2)',[PopSize,1])+new_POPS(1:PopSize,1:SearchDimension))./2);
end

function [new_POPS, new_POPS_indices, B] = selectAndUpdatePopulation_uniquePOPM(POPM, B, PopSize, SearchDimension, selection_prob)
    % Get unique POPM
    [uniquePOPM, ~, ic] = unique(POPM(:,1:SearchDimension+1), 'rows', 'stable');
    
    % Initialize
    valid_indices = true(size(uniquePOPM, 1), 1);  % For marking valid individuals
    new_POPS = zeros(PopSize, size(POPM, 2));
    new_POPS_indices = zeros(PopSize, 1);
    new_POPS_count = 0;

    while new_POPS_count < PopSize && any(valid_indices)
        % Only select from valid individuals
        valid_selection_prob = selection_prob .* valid_indices;
        valid_selection_prob = valid_selection_prob / sum(valid_selection_prob);
        
        % Use randsample function to select individuals
        selected_index = randsample(size(uniquePOPM, 1), 1, true, valid_selection_prob);
        selected_individual = uniquePOPM(selected_index, :);
        
        % Check the number of current selected individuals in the new population
        duplicate_count = sum(all(new_POPS(1:new_POPS_count, 1:SearchDimension) == selected_individual(1:SearchDimension), 2));
        
        % If count doesn't exceed the limit, add to new population
        if duplicate_count < max_duplicates
            new_POPS_count = new_POPS_count + 1;
            new_POPS(new_POPS_count, :) = selected_individual;
            % Find index in original POPM
            original_indices = find(ic == selected_index);
            new_POPS_indices(new_POPS_count) = original_indices(randi(length(original_indices)));
        else
            % Mark all individuals identical to the selected one as invalid
            valid_indices(selected_index) = false;
        end
    end

    % If new population is not full, fill remaining positions with valid individuals from original population
    if new_POPS_count < PopSize
        remaining_indices = find(valid_indices);
        remaining_count = PopSize - new_POPS_count;
        if isempty(remaining_indices)
            % If not enough valid individuals, reset all individuals as valid
            remaining_indices = 1:size(uniquePOPM, 1);
        end
        fill_indices = remaining_indices(randperm(length(remaining_indices), remaining_count));
        for i = 1:length(fill_indices)
            new_POPS_count = new_POPS_count + 1;
            new_POPS(new_POPS_count, :) = uniquePOPM(fill_indices(i), :);
            original_indices = find(ic == fill_indices(i));
            new_POPS_indices(new_POPS_count) = original_indices(randi(length(original_indices)));
        end
    end

    % Ensure best individual is in new_POPS_indices
    % Find index of individual with lowest fitness in POPM
    [~, best_index] = min(POPM(:, SearchDimension+1));
    
    % Check if best individual is in new_POPS_indices
    if ~ismember(best_index, new_POPS_indices)
        % Calculate occurrence count of each index in new_POPS_indices
        [unique_indices, ~, ic] = unique(new_POPS_indices);
        index_counts = accumarray(ic, 1);
        
        % Find all indices with maximum occurrence count
        max_count = max(index_counts);
        most_frequent_indices = unique_indices(index_counts == max_count);
        
        % Randomly select one from these indices
        most_frequent_index = most_frequent_indices(randi(length(most_frequent_indices)));
        
        % Find positions of most frequently occurring index in new_POPS_indices
        replace_positions = find(new_POPS_indices == most_frequent_index);
        
        % Randomly select a position to replace
        replace_position = replace_positions(randi(length(replace_positions)));
        
        % Replace the selected position with best individual's index
        new_POPS_indices(replace_position) = best_index;
        
        % Update corresponding individual in new_POPS
        new_POPS(replace_position, :) = POPM(best_index, :);
    end

    % Return indices of unselected individuals in POPM
    unselected_indices = setdiff(1:size(POPM, 1), new_POPS_indices);

    % Add unselected individuals to archive B
    B = updateArchive(B, POPM, unselected_indices, SearchDimension, PopSize);
end

end