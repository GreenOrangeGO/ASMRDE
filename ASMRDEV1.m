function [Result,AdaptFuncValue] = ASMRDEV1(PopSize,SearchDimension,SScope,FuncNo,MaxGen)

% 定义随机数种子
rand('state',sum(100*clock));
% 针对CEC测试PopSize_Input等于问题维度
PopSize_Input = SearchDimension;
% 定义最大函数评价次数
NFEmax = MaxGen*PopSize_Input;
% 定义初始种群
POPS = rand(PopSize,SearchDimension+1);
% 随机初始化种群
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
% 针对CEC测试问题，FuncNo为函数序号
POPS(1:PopSize,SearchDimension+1) = cec17_func(POPS(1:PopSize,1:SearchDimension)',FuncNo)';
% 找到带有最优适应度值的个体
[~,row] = min(POPS(1:PopSize,SearchDimension+1));
% 定义BEST为最有个体
BEST = POPS(row,1:SearchDimension+1);
% 记录被评价的个体的适应度值，以便后续绘制迭代曲线
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
% 初始时刻评级次数
nfes = PopSize;

% 定义初始时刻迭代次数
k = 1;

% 参数
CR = 0.7;

T = 9000/k;

Stagnant = SearchDimension*0.3;

%初始化存档B
B = [];

% 参数F的变化范围
delta = 0.1;

% 在主循环开始前定义一个函数来计算多样性
function diversity = calculateDiversity(population)
    [uniquePopulation, ~, ic] = unique(population, 'rows');
    diversity = size(uniquePopulation, 1);
end

function selection_prob = calculateSelectionProbability(POPS, PopSize, SearchDimension, T)
    % 计算适应度（直接使用原始适应度值，较小的值更好）
    fitness = -POPS(1:PopSize, SearchDimension+1);

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
        exp_fitness = exp(fitness_normalized / T);
        selection_prob = exp_fitness / sum(exp_fitness);
        
        % 添加小的正数 ε 以确保没有零概率
        epsilon = 1e-10;
        selection_prob = selection_prob + epsilon;
        selection_prob = selection_prob / sum(selection_prob);
    end
end

while nfes < NFEmax
    
    % 更新迭代次数
    k = k + 1;

    %显示当前最佳适应度和多样性 (每5000次迭代输出一次)
    if mod(k, 5000) == 0 || k == 1 || k == MaxGen
        BestFitness = min(POPS(1:PopSize,SearchDimension+1));
        POPSdiversity = calculateDiversity(POPS(1:PopSize,1:SearchDimension));
        Bdiversity = calculateDiversity(B(1:PopSize,1:SearchDimension));
        POPMdiversity = calculateDiversity(POPM(1:PopSize*2,1:SearchDimension));
        str=['Function is F' num2str(FuncNo) ' Iterations:' num2str(k)...
            ' Best fitness:' num2str(BestFitness-FuncNo*100,'%e')...
            ' Population diversity:' num2str(POPSdiversity)...
            ' B diversity:' num2str(Bdiversity)...
            ' POPM diversity:' num2str(POPMdiversity)];
        disp(str);
    end

    %500代选择优良个体一次
    if mod(k,500) ==0 

        % 初始化停滞计数器
        stagnantCounter = zeros(PopSize, 1);
        lastGeneration = zeros(PopSize, 1);

        % 4.3 选择操作
        selection_prob = calculateSelectionProbability(POPS, PopSize, SearchDimension, T);

        % 设置最大相同个体数量
        max_duplicates = floor(PopSize / 3);

        % 初始化
        valid_indices = true(PopSize, 1);  % 用于标记有效的个体
        new_POPS = zeros(size(POPS));
        new_POPS_indices = zeros(PopSize, 1);
        new_POPS_count = 0;

        while new_POPS_count < PopSize && any(valid_indices)
            % 只从有效的个体中选择
            valid_selection_prob = selection_prob .* valid_indices;
            valid_selection_prob = valid_selection_prob / sum(valid_selection_prob);
            
            % 使用 randsample 函数选择个体
            selected_index = randsample(PopSize, 1, true, valid_selection_prob);
            selected_individual = POPS(selected_index, :);
            
            % 检查当前选中个体在新种群中的数量
            duplicate_count = sum(all(new_POPS(1:new_POPS_count, 1:SearchDimension) == selected_individual(1:SearchDimension), 2));
            
            % 如果数量未超过限制，则添加到新种群
            if duplicate_count < max_duplicates
                new_POPS_count = new_POPS_count + 1;
                new_POPS(new_POPS_count, :) = selected_individual;
                new_POPS_indices(new_POPS_count) = selected_index;
            else
                % 标记所有与选中个体相同的个体为无效
                duplicate_indices = find(all(POPS(:, 1:SearchDimension) == selected_individual(1:SearchDimension), 2));
                valid_indices(duplicate_indices) = false;
            end
        end

        % 如果新种群没有填满，用原始种群中的有效个体填充剩余位置
        if new_POPS_count < PopSize
            remaining_indices = find(valid_indices);
            remaining_count = PopSize - new_POPS_count;
            if length(remaining_indices) == 0
                % 如果有效个体不足，重置所有个体为有效
                remaining_indices = 1:PopSize;
            end
            fill_indices = remaining_indices(randperm(length(remaining_indices), remaining_count));
            new_POPS(new_POPS_count+1:PopSize, :) = POPS(fill_indices, :);
            new_POPS_indices(new_POPS_count+1:PopSize) = fill_indices;
        end

        % 返回solution中未被选中个体的index
        unselected_indices = setdiff(1:PopSize, new_POPS_indices);
        
        % 将未被选中的个体添加到存档B中
        B = [B; POPS(unselected_indices, :)];

        % 如果B的大小超过PopSize，保留最好的PopSize个个体
        if size(B, 1) > PopSize
            [~, sorted_indices] = sort(B(:, end));
            B = B(sorted_indices(1:PopSize), :);
        end

        % 建立种群并集P∪A
        POPM = [POPS;B];

        % 获得并集个体数
        PopMSolution = size(POPM,1);

        % % 输出命令行
        % [BestFitness, ~] = min(POPM(1:end,SearchDimension+1));
        % str = ['F' num2str(FuncNo) ' Iter:' num2str(k)...
        %     ' Eval:' num2str(nfes) ' Best:'...
        %     num2str(BestFitness-100*FuncNo,'%e')];
        % disp(str);

        % 挑选两个随机数并且不等于当前个体
        nrandI = 2;
        rr = zeros(PopSize, nrandI);
        for i = 1:PopSize
            original_index = new_POPS_indices(i);
            available_index = [1:original_index-1 original_index+1:size(POPM,1)];
            sequence_index = randperm (numel(available_index));
            rr(i,:) = available_index(sequence_index(1:nrandI));
        end

        % 计算F取值
        MeanF = 0.9 - 0.8 * (nfes-1)/(NFEmax-1);
        F = unifrnd(MeanF - delta,MeanF + delta);

        % 选择Xpbest基向量，首先计算参数p：精英比率，使其按照线性递减变化
        [~,indexbest] = sort(POPM(1:end,1+SearchDimension),'ascend');
        p = 1 - 1*(nfes/NFEmax)^1;
        pNP = round(max(1,PopMSolution*p));

        % 产生随机数挑选Xpbest
        randindex = ceil(rand(PopSize,1)*pNP);
        randindex = indexbest(randindex);
        Xpbest = POPM(randindex,1:SearchDimension);

        % 初始化变异向量V
        V = zeros(PopSize,SearchDimension);

        % 利用公式DE/current-to-pbest/1 with external archive B产生V
        V(1:PopSize,1:SearchDimension) = new_POPS(1:PopSize,1:SearchDimension) + F * (Xpbest - new_POPS(1:PopSize,1:SearchDimension)) +...
            F*(POPM(rr(1:PopSize,1),1:SearchDimension) - POPM(rr(1:PopSize,2),1:SearchDimension));

        % 将V中超出边界的维度进行处理，处理方式为取对应X和边界范围的平均数
        V(1:PopSize,1:SearchDimension) = ...
            ((V(1:PopSize,1:SearchDimension)>=repmat(SScope(:,1)',[PopSize,1]))&(V(1:PopSize,1:SearchDimension)<=repmat(SScope(:,2)',[PopSize,1]))).*(V(1:PopSize,1:SearchDimension))+...
            (V(1:PopSize,1:SearchDimension)<repmat(SScope(:,1)',[PopSize,1])).*((repmat(SScope(:,1)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2)+...
            (V(1:PopSize,1:SearchDimension)>repmat(SScope(:,2)',[PopSize,1])).*((repmat(SScope(:,2)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2);
        
        POPS = new_POPS;
    else
        if k > 500
            % 检查stagnantCounter中是否有元素等于Stagnant
            replace_indices = find(stagnantCounter == Stagnant);
    
            if ~isempty(replace_indices)

                % 在替换循环开始前，保存原始的停滞个体
                original_stagnant_individuals = POPS(replace_indices, :);
                
                selection_prob = calculateSelectionProbability(B, size(B,1), SearchDimension, T);
                
                % 设置最大相同个体数量
                max_duplicates = floor(PopSize / 4);
    
                % 初始化替换个体数组和有效索引
                replacement_individuals = zeros(length(replace_indices), size(POPS, 2));
                replaced_count = 0;
                valid_indices = true(size(B, 1), 1);  % 用于标记B中有效的个体
    
                while replaced_count < length(replace_indices) && any(valid_indices)
                    % 只从有效的个体中选择
                    valid_selection_prob = selection_prob .* valid_indices;
                    valid_selection_prob = valid_selection_prob / sum(valid_selection_prob);
                    
                    % 从B中选择一个个体
                    selected_index = randsample(size(B, 1), 1, true, valid_selection_prob);
                    selected_individual = B(selected_index, :);
                    
                    % 检查POPS中与选中个体相同的个体数量
                    duplicate_count = sum(all(POPS(:, 1:SearchDimension) == selected_individual(1:SearchDimension), 2));
                    
                    % 如果数量未超过限制，则添加到替换个体数组
                    if duplicate_count < max_duplicates
                        replaced_count = replaced_count + 1;
                        replacement_individuals(replaced_count, :) = selected_individual;
                        
                        % 更新POPS以反映新添加的个体（这步是为了准确计算后续的duplicate_count）
                        POPS(replace_indices(replaced_count), :) = selected_individual;
                    else
                        % 标记所有与选中个体相同的个体为无效
                        duplicate_indices = find(all(B(:, 1:SearchDimension) == selected_individual(1:SearchDimension), 2));
                        valid_indices(duplicate_indices) = false;
                    end
                end
    
                % 如果替换个体不足，用B中的其他有效个体填充
                if replaced_count < length(replace_indices)
                    remaining_indices = find(valid_indices);
                    remaining_count = length(replace_indices) - replaced_count;
                    if length(remaining_indices) == 0
                        % 如果有效个体不足，重置所有个体为有效
                        remaining_indices = 1:size(B, 1);
                    end
                    fill_indices = remaining_indices(randperm(length(remaining_indices), remaining_count));
                    replacement_individuals(replaced_count+1:end, :) = B(fill_indices, :);
                end
    
                % 非常不好，全面退化，停滞个体似乎没有一根毛的价值（观察
    
                % 将原始的停滞个体添加到存档B中
                B = [B; original_stagnant_individuals];
                
                % 如果B的大小超过PopSize，保留最好的PopSize个个体
                if size(B, 1) > PopSize
                    [~, sorted_indices] = sort(B(:, end));
                    B = B(sorted_indices(1:PopSize), :);
                end
    
                % 替换POPS中对应的个体
                POPS(replace_indices, :) = replacement_individuals;
    
                % 重置被替换个体的stagnantCounter
                stagnantCounter(replace_indices) = 0;
            end
        end
        %正常循环
        % 建立种群并集P∪B
        POPM = [POPS;B];

        % 获得并集个体数
        PopMSolution = size(POPM,1);

        % % 输出命令行
        % [BestFitness, ~] = min(POPM(1:end,SearchDimension+1));
        % str = ['F' num2str(FuncNo) ' Iter:' num2str(k)...
        %     ' Eval:' num2str(nfes) ' Best:'...
        %     num2str(BestFitness-100*FuncNo,'%e')];
        % disp(str);

        % 挑选两个随机数并且不等于当前个体
        nrandI = 2;
        rr = zeros(PopSize, nrandI);
        for i = 1:PopSize
            available_index = [1:i-1 i+1:size(POPM,1)];
            sequence_index = randperm (numel(available_index));
            rr(i,:) = available_index(sequence_index(1:nrandI));
        end

        % 计算F取值
        MeanF = 0.9 - 0.8 * (nfes-1)/(NFEmax-1);
        F = unifrnd(MeanF - delta,MeanF + delta);

        % 选择Xpbest基向量，首先计算参数p：精英比率，使其按照线性递减变化
        [~,indexbest] = sort(POPM(1:end,1+SearchDimension),'ascend');
        p = 1 - 1*(nfes/NFEmax)^1;
        pNP = round(max(1,PopMSolution*p));

        % 产生随机数挑选Xpbest
        randindex = ceil(rand(PopSize,1)*pNP);
        randindex = indexbest(randindex);
        Xpbest = POPM(randindex,1:SearchDimension);

        % 初始化变异向量V
        V = zeros(PopSize,SearchDimension);

        % 利用公式DE/current-to-pbest/1 with external archive B产生V
        V(1:PopSize,1:SearchDimension) = POPS(1:PopSize,1:SearchDimension) + F * (Xpbest - POPS(1:PopSize,1:SearchDimension)) +...
            F*(POPM(rr(1:PopSize,1),1:SearchDimension) - POPM(rr(1:PopSize,2),1:SearchDimension));

        % 将V中超出边界的维度进行处理，处理方式为取对应X和边界范围的平均数
        V(1:PopSize,1:SearchDimension) = ...
            ((V(1:PopSize,1:SearchDimension)>=repmat(SScope(:,1)',[PopSize,1]))&(V(1:PopSize,1:SearchDimension)<=repmat(SScope(:,2)',[PopSize,1]))).*(V(1:PopSize,1:SearchDimension))+...
            (V(1:PopSize,1:SearchDimension)<repmat(SScope(:,1)',[PopSize,1])).*((repmat(SScope(:,1)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2)+...
            (V(1:PopSize,1:SearchDimension)>repmat(SScope(:,2)',[PopSize,1])).*((repmat(SScope(:,2)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2);
    
    end
    % 定义试验向量U
    U = POPS(1:PopSize,1:SearchDimension+1);

    % 随机选择一维度，这一维必定交叉
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);

    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);

    % I为逻辑矩阵
    I = (rand(PopSize,SearchDimension) < repmat(CR,[PopSize,SearchDimension])) | (j == jRand);
    U(I) = V(I);

    % 评价试验向量
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',FuncNo);

    % 更新评价次数
    nfes = nfes + PopSize;

    tmp = (U(1:PopSize,SearchDimension+1) < POPS(1:PopSize,SearchDimension+1));        % 产生 tmp 作为子代替换父代的标记
        
    if k > 500
        % Update stagnantCounter based on conditions
        stagnantCounter = stagnantCounter + ((lastGeneration == 0 & tmp == 0) | (lastGeneration == 1 & tmp == 0));
        stagnantCounter(lastGeneration == 0 & tmp == 1) = 0;
        lastGeneration = tmp;
    end

    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp 的D维复制版本

    % 更新种群
    POPS(1:PopSize,1:SearchDimension+1) = temp.*U(1:PopSize,1:SearchDimension+1) + (1-temp).*POPS(1:PopSize,1:SearchDimension+1);

    % 记录函数评价值
    for i = 1:PopSize
        SearchProcess(1,nn) = ...
            (SearchProcess(1,nn-1) > POPS(i,SearchDimension+1)).*POPS(i,SearchDimension+1) +...
            (SearchProcess(1,nn-1) <= POPS(i,SearchDimension+1)).*SearchProcess(1,nn-1);
        nn = nn + 1;
    end

    % 更新P∪U
    POPM = [POPS;B];

    % 更新合并种群
    POPM = [POPS(1:end,:);B];

    % 从合并种群中挑出最好的个体并记录
    [~,row] = min(POPM(1:end,SearchDimension+1));
    BEST = POPM(row,1:SearchDimension+1);

end
% 定义绘制迭代曲线代数间隔
kk = 1:MaxGen;

% 输出结果AdaptFuncValue用于绘制曲线
AdaptFuncValue = SearchProcess(PopSize_Input.*kk);

% 输出结果Result用于计算Mean、SD、Best、Worst值，随后用于统计分析
Result = BEST;
end
