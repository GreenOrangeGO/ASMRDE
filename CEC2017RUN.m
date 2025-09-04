clear;
clc;

% Get the number of available logical cores
numCores = feature('numcores');
% numCores = 13;

% Enable parallel computing pool, try to use all available logical cores
if isempty(gcp('nocreate'))
    parpool('local', numCores);
end

% Display the actual number of workers used
p = gcp;
fprintf('Using %d workers in the parallel pool\n', p.NumWorkers);

% Parameters
SearchDimension = 10;

if SearchDimension == 10
    PopSize = 90;
elseif SearchDimension == 30
    PopSize = 150;
elseif SearchDimension == 50
    PopSize = 210;
elseif SearchDimension == 100
    PopSize = 260;
end

SearchScope = repmat([-100 100], SearchDimension, 1);
LoopCount = ceil(SearchDimension*10000/PopSize);
RunCount = 51;

% Define the list of algorithms to run
algorithms = {'ASMRDE'}; % More algorithms can be added
current_algorithm = 1;

while current_algorithm <= length(algorithms)
    alg_name = algorithms{current_algorithm};

    % Check if results file already exists
    result_file = [alg_name '_17D' num2str(SearchDimension) '.mat'];
    if exist(result_file, 'file')
        fprintf('Results for %s already exist, skipping...\n', alg_name);
        current_algorithm = current_algorithm + 1;
        continue;
    end

    fprintf('Running algorithm: %s\n', alg_name);

    for FuncNum = [1,3:30]

        Result_temp = zeros(RunCount, SearchDimension+1);
        AdaptFuncValue_temp = zeros(RunCount, 10000);

        tic
        parfor pp = 1:RunCount
            [Result_temp(pp,:), AdaptFuncValue_temp(pp,:)] = feval(alg_name, PopSize, SearchDimension, SearchScope, FuncNum, LoopCount);
        end
        toc

        % Calculate statistics
        % Create variables named by function number to store data
        eval(['AdaptFuncValue_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = AdaptFuncValue_temp;']);
        eval(['AdaptFuncValue_avg_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = mean(AdaptFuncValue_temp, 1);']);

        eval(['Result_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = Result_temp;']);
        eval(['Result_avg_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = mean(Result_temp, 1);']);

        % Calculate and display statistical results
        format shortE; eval(['Mean_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = mean(Result_temp(:,SearchDimension+1)) - ' num2str(FuncNum*100) ';']);
        format shortE; eval(['Best_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = min(Result_temp(:,SearchDimension+1)) - ' num2str(FuncNum*100) ';']);
        format shortE; eval(['Worst_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = max(Result_temp(:,SearchDimension+1)) - ' num2str(FuncNum*100) ';']);
        format shortE; eval(['SD_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name ' = std(Result_temp(:,SearchDimension+1) - ' num2str(FuncNum*100) ', 1);']);

        % Output calculation results
        fprintf('Function name: f%d\n', FuncNum);
        fprintf('Mean: %.6e\n', eval(['Mean_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name]));
        fprintf('Best: %.6e\n', eval(['Best_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name]));
        fprintf('Worst: %.6e\n', eval(['Worst_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name]));
        fprintf('Standard Deviation: %.6e\n', eval(['SD_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name]));
        fprintf('\n');
    end

    % Save results
    filename = [alg_name '_17D' num2str(SearchDimension) '.mat'];
    save(filename);

    % Collect mean values and write to Excel
    funcNums = [1,3:30];
    MeanValues = zeros(length(funcNums), 1);

    for i = 1:length(funcNums)
        FuncNum = funcNums(i);
        MeanValues(i) = eval(['Mean_f' num2str(FuncNum) '_D' num2str(SearchDimension) '_' alg_name]);
    end

    % Write to Excel
    writematrix(MeanValues, [alg_name '_17D' num2str(SearchDimension) '.xlsx']);
    
    % Clean up current algorithm variables, but keep important program control variables
    vars_to_keep = {'current_algorithm', 'algorithms', 'PopSize', 'SearchDimension', ...
                    'SearchScope', 'LoopCount', 'RunCount', 'p', 'numCores'};

    % Get all workspace variables
    all_vars = who;

    % Find variables to clear
    vars_to_clear = setdiff(all_vars, vars_to_keep);

    % Clear variables
    clear(vars_to_clear{:});

    % Continue to next algorithm
    current_algorithm = current_algorithm + 1;
end
