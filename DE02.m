function [Result,AdaptFuncValue] = DE02(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)

% 1. 输入参数验证
if nargin < 4
    error('输入参数的个数错误。')
end

[row,colum] = size(PopSize);
if row>1 || colum>1
    error('输入的种群大小错误，应该是一个 1 行 1 列的整数。');
end
if PopSize ~=fix (PopSize)
    error('输入的种群大小错误，应该是一个整数。');    
end
if PopSize < 6
    error('输入的种群大小错误，应该是一个大于或等于 6 的整数。');    
end

[row,colum] = size(SearchDimension);
if row>1 || colum>1
    error('输入的搜索空间维数错误，应该是一个 1 行 1 列的整数。');
end
if SearchDimension ~= fix(SearchDimension)
    error('输入的搜索空间维数错误，应该是一个整数。');    
end
if SearchDimension < 1
    error('输入的搜索空间维数错误，应该是一个大于或等于 1 的整数。');    
end

[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('输入的各维取值范围错误，应该是一个 %g 行 2 列的实数矩阵。',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('输入的各维取值范围错误，第 %g 行中的第二个数据应该大于第一个数据。',d);
    end
end

if nargin == 4               
    LoopCount = 1000;  
end

[row,colum] = size(LoopCount);
if row>1 || colum>1
    error('输入的总迭代次数错误，应该是一个 1 行 1 列的整数。');
end
if LoopCount ~= fix(LoopCount)
    error('输入的总迭代次数错误，应该是一个整数。');
end
if LoopCount < 1
    error('输入的总迭代次数错误，应该是一个大于或等于 1 的整数。');
end   

if nargout ~= 2
    error('输出参数的个数错误。')
end

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

% 4. 主循环
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
    % 4.1 显示当前最佳适应度 
    BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Best fitness:'...
        num2str(BestFitness-AdaptFunc*100,'%e')];
    disp(str);
    nfes = nfes + PopSize;
    
    % 4.2 设置 DE 参数
    F = 0.5;
    CR = 0.9;

    % 4.3 选择操作
    % 计算适应度（取负值，使得原本较小的值变大）
    fitness = -Solution(1:PopSize, SearchDimension+1);

    % 标准化适应度（Z-score标准化）
    fitness_mean = mean(fitness);
    fitness_std = std(fitness);
    fitness_normalized = (fitness - fitness_mean) / fitness_std;

    % 设置温度参数
    T = 1;

    % 使用softmax计算选择概率
    exp_fitness = exp(fitness_normalized / T);
    selection_prob = exp_fitness / sum(exp_fitness);

    % 基于选择概率随机选择新种群
    new_population_indices = randsample(PopSize, PopSize, true, selection_prob);
    new_population = Solution(new_population_indices, :);

    % 4.4 生成随机索引
    nrandI = 5;
    r = zeros(PopSize, nrandI);
    for i = 1:PopSize
        original_index = new_population_indices(i);
        available_indices = setdiff(1:PopSize, original_index);
        r(i,:) = available_indices(randperm(PopSize-1, nrandI));
    end

    % 4.5 变异操作
    V(1:PopSize,1:SearchDimension) = new_population(:,1:SearchDimension) + F * (Solution(r(:,2),1:SearchDimension) - Solution(r(:,3),1:SearchDimension));

    % 4.6 边界处理
    V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + new_population(:,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + new_population(:,1:SearchDimension))./2);

    % 4.7 交叉操作
    U = new_population;
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(:,1:SearchDimension) = I .* V(:,1:SearchDimension) + (~I) .* new_population(:,1:SearchDimension);

    % 4.8 评估新解
    U(:,SearchDimension+1) = cec17_func(U(:,1:SearchDimension)',AdaptFunc); 

    % 4.9 选择
    tmp = (U(:,SearchDimension+1) < new_population(:,SearchDimension+1));
    temp = repmat(tmp,1,SearchDimension+1);
    Solution(1:PopSize,:) = temp.*U + (1-temp).*new_population;

    % 4.10 更新最佳解
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);

    % 4.11 记录当前最佳适应度值
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
end

% 5. 返回结果
Result = Solution(PopSize+1,:);

end



