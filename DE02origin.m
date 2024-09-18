function [Result,AdaptFuncValue] = DE02(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)

% 2. 初始化
rand('state',sum(100*clock)); 

Solution = rand(PopSize+1,SearchDimension+1);
A = generator1(PopSize,SearchDimension);
Solution(1:PopSize,1:SearchDimension) = A;
for d = 1:SearchDimension
    if SearchScope(d,1)==-inf && SearchScope(d,2)~=inf
        Solution(:,d) = unifrnd(SearchScope(d,2)-10^10,SearchScope(d,2),PopSize,1);
    end
    
    if SearchScope(d,1)~=-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(SearchScope(d,1),SearchScope(d,1)+10^10,PopSize,1);
    end
    
    if SearchScope(d,1)==-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(-10^10,10^10,PopSize,1);
    end
    
    if SearchScope(d,1)~=-inf && SearchScope(d,2)~=inf
        Solution(:,d) = Solution(:,d)*(SearchScope(d,2)-SearchScope(d,1))+SearchScope(d,1);
    end                                                                                                                                               
end

% 3. 评估初始种群
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc); 

[~,row] = min(Solution(1:PopSize,SearchDimension+1));  
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);

% 在主循环开始前定义一个函数来计算多样性
function diversity = calculateDiversity(population)
    [uniquePopulation, ~, ic] = unique(population, 'rows');
    diversity = size(uniquePopulation, 1);
end

% 4. 主循环
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
    % 4.1 显示当前最佳适应度和多样性 (每5000次迭代输出一次)
    if mod(k, 5000) == 0 || k == 1 || k == LoopCount
        BestFitness = min(Solution(1:PopSize,SearchDimension+1));
        diversity = calculateDiversity(Solution(1:PopSize,1:SearchDimension));
        str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
            ' Best fitness:' num2str(BestFitness-AdaptFunc*100,'%e')...
            ' Population diversity:' num2str(diversity)];
        disp(str);
    end
    nfes = nfes + PopSize;
    
    % 4.2 设置 DE 参数
    F = 0.5;
    CR = 0.9;

    % 设置温度参数
    T = 60.1-0.006*k;

    %500代选择优良个体一次
    if mod(k,500) ==0 
    
        %初始化存档B
        B = [];

        stagnantCounter = zeros(PopSize, 1);

        % 4.3 选择操作
        % 计算适应度（直接使用原始适应度值，较小的值更好）
        fitness = -Solution(1:PopSize, SearchDimension+1);

        % 标准化适应度（Z-score标准化），但要处理标准差为0的情况
        fitness_mean = mean(fitness);
        fitness_std = std(fitness);

        if fitness_std == 0
            % 如果标准差为0，说明所有适应度值相同
            % 在这种情况下，我们可以给每个个体相等的选择概率
            selection_prob = ones(PopSize, 1) / PopSize;
        else
            % 正常进行标准化和后续计算
            fitness_normalized = (fitness - fitness_mean) / fitness_std;
            
            % 使用softmax计算选择概率
            exp_fitness = exp(fitness_normalized/ T);
            selection_prob = exp_fitness / sum(exp_fitness);
            
            % 添加小的正数 ε 以确保没有零概率
            epsilon = 1e-10;
            selection_prob = selection_prob + epsilon;
            selection_prob = selection_prob / sum(selection_prob);
        end

        % 使用 randsample 函数
        new_population_indices = randsample(PopSize, PopSize, true, selection_prob);
        new_population = Solution(new_population_indices, :);

        % 返回solution中未被选中个体的index
        unselected_indices = setdiff(1:PopSize, new_population_indices);
        
        % 将未被选中的个体添加到存档B中
        B = [B; Solution(unselected_indices, :)];

        % 如果B的大小超过PopSize，保留最好的PopSize个个体
        if size(B, 1) > PopSize
            [~, sorted_indices] = sort(B(:, end));
            B = B(sorted_indices(1:PopSize), :);
        end

        % 4.4 生成随机索引
        nrandI = 5;
        r = zeros(PopSize, nrandI);
        for i = 1:PopSize
            original_index = new_population_indices(i);
            available_indices = setdiff(1:PopSize, original_index);
            r(i,:) = available_indices(randperm(PopSize-1, nrandI));
        end

        % 4.5 变异操作
        % 优化：使用矩阵运算代替循环
        V = Solution(r(1:PopSize,1),1:SearchDimension) + F * (Solution(r(:,2),1:SearchDimension) - Solution(r(:,3),1:SearchDimension));

        % 4.6 边界处理
        % 优化：使用逻辑索引代替条件语句
        lower_bound = repmat(SearchScope(:,1)', [PopSize,1]);
        upper_bound = repmat(SearchScope(:,2)', [PopSize,1]);
        V_pop = new_population(:,1:SearchDimension);

        V_out_lower = V < lower_bound;
        V_out_upper = V > upper_bound;
        V_in_bounds = ~(V_out_lower | V_out_upper);

        V = V .* V_in_bounds + ...
            ((lower_bound + V_pop) / 2) .* V_out_lower + ...
            ((upper_bound + V_pop) / 2) .* V_out_upper;

        % 4.7 交叉操作
        % 优化：使用矩阵运算代替循环
        jRand = repmat(ceil(rand(PopSize,1)*SearchDimension), 1, SearchDimension);
        j = repmat(1:SearchDimension, PopSize, 1);
        I = (rand(PopSize,SearchDimension) < CR) | (j == jRand);
        U(:,1:SearchDimension) = I .* V + (~I) .* new_population(:,1:SearchDimension);

        % 4.8 评估新解
        U(:,SearchDimension+1) = cec17_func(U(:,1:SearchDimension)',AdaptFunc); 

        % 4.9 选择
        tmp = (U(:,SearchDimension+1) < new_population(:,SearchDimension+1));
        temp = repmat(tmp,1,SearchDimension+1);
        Solution(1:PopSize,:) = temp.*U + (1-temp).*new_population; 

    else
        
        if k > 500
        % 检查stagnantCounter中是否有元素大于32
        replace_indices = find(stagnantCounter > 32);

            if ~isempty(replace_indices)
            % 从存档B中随机选择个体
            % 选择操作
            % 计算适应度（直接使用原始适应度值，较小的值更好）
            fitness = -B(:, SearchDimension+1);

            % 标准化适应度（Z-score标准化），但要处理标准差为0的情况
            fitness_mean = mean(fitness);
            fitness_std = std(fitness);

                if fitness_std == 0
                    % 如果标准差为0，说明所有适应度值相同
                    % 在这种情况下，我们可以给每个个体相等的选择概率
                    selection_prob = ones(height(B), 1) / height(B);
                else
                    % 正常进行标准化和后续计算
                    fitness_normalized = (fitness - fitness_mean) / fitness_std;

                    % 使用softmax计算选择概率
                    exp_fitness = exp(fitness_normalized/ T);
                    selection_prob = exp_fitness / sum(exp_fitness);

                    % 添加小的正数 ε 以确保没有零概率
                    epsilon = 1e-10;
                    selection_prob = selection_prob + epsilon;
                    selection_prob = selection_prob / sum(selection_prob);
                end

                % 使用 randsample 函数
                replace_population_indices = randsample(height(B), length(replace_indices), true, selection_prob);
                replacement_individuals = B(replace_population_indices, :);

                %非常不好，全面退化，停滞个体似乎没有一根毛的价值

                % % 将被替换的个体添加到存档B中
                % B = [B; Solution(replace_indices, :)];
                % 
                % % 如果B的大小超过PopSize，保留最好的PopSize个个体
                % if size(B, 1) > PopSize
                %     [~, sorted_indices] = sort(B(:, end));
                %     B = B(sorted_indices(1:PopSize), :);
                % end

                % 替换Solution中对应的个体
                Solution(replace_indices, :) = replacement_individuals;

                % 重置被替换个体的stagnantCounter
                stagnantCounter(replace_indices) = 0;
            end
        end

        %正常循环
        nrandI = 5;
        r = zeros(PopSize,nrandI);
        for i = 1:PopSize
            available_index = [1:i-1 i+1:PopSize];
            sequence_index = randperm (numel(available_index));
            r(i,:) = available_index(sequence_index(1:nrandI));
        end
    
        V(1:PopSize,1:SearchDimension) = Solution(r(1:PopSize,1),1:SearchDimension) + F *(Solution(r(1:PopSize,2),1:SearchDimension) - Solution(r(1:PopSize,3),1:SearchDimension));
    
        % 边界修正
        V(1:PopSize,1:SearchDimension) = ...
         ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
         (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
         (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

        % 以下为基于矩阵运算的交叉操作
        U = Solution(1:PopSize,:);
        jRand = ceil(rand(PopSize,1)*SearchDimension);
        jRand = repmat(jRand,[1,SearchDimension]);
        j = 1:SearchDimension;
        j = repmat(j,[PopSize,1]);
        I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
        U(I) = V(I);
        
        U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc); 
        tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % 产生 tmp 作为子代替换父代的标记
        
        if k > 500
            lastGeneration = tmp;
            % Update stagnantCounter based on conditions
            stagnantCounter = stagnantCounter + ((lastGeneration == 0 & tmp == 0) | (lastGeneration == 1 & tmp == 0));
            stagnantCounter(lastGeneration == 0 & tmp == 1) = 0;
        end

        temp = repmat(tmp,1,SearchDimension+1);                                                % tmp 的D维复制版本
        Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % 子代成功的，保留U；子代不成功的，保留父代
    end

    % % 4.1 显示当前最佳适应度和多样性 (每5000次迭代输出一次)
    % if mod(k, 5000) == 0 || k == 1 || k == LoopCount
    %     BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    %     diversity = calculateDiversity(Solution(1:PopSize,1:SearchDimension));
    %     str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
    %         ' Best fitness:' num2str(BestFitness-AdaptFunc*100,'%e')...
    %         ' Population diversity:' num2str(diversity)];
    %     disp(str);
    % end
    
    % 4.10 更新最佳解
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);

    % 4.11 记录当前最佳适应度值
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
end

% 5. 返回结果
Result = Solution(PopSize+1,:);

end