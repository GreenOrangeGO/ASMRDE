function [Result,AdaptFuncValue] = DE02(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)

% 1. ���������֤
if nargin < 4
    error('��������ĸ�������')
end

[row,colum] = size(PopSize);
if row>1 || colum>1
    error('�������Ⱥ��С����Ӧ����һ�� 1 �� 1 �е�������');
end
if PopSize ~=fix (PopSize)
    error('�������Ⱥ��С����Ӧ����һ��������');    
end
if PopSize < 6
    error('�������Ⱥ��С����Ӧ����һ�����ڻ���� 6 ��������');    
end

[row,colum] = size(SearchDimension);
if row>1 || colum>1
    error('����������ռ�ά������Ӧ����һ�� 1 �� 1 �е�������');
end
if SearchDimension ~= fix(SearchDimension)
    error('����������ռ�ά������Ӧ����һ��������');    
end
if SearchDimension < 1
    error('����������ռ�ά������Ӧ����һ�����ڻ���� 1 ��������');    
end

[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('����ĸ�άȡֵ��Χ����Ӧ����һ�� %g �� 2 �е�ʵ������',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('����ĸ�άȡֵ��Χ���󣬵� %g ���еĵڶ�������Ӧ�ô��ڵ�һ�����ݡ�',d);
    end
end

if nargin == 4               
    LoopCount = 1000;  
end

[row,colum] = size(LoopCount);
if row>1 || colum>1
    error('������ܵ�����������Ӧ����һ�� 1 �� 1 �е�������');
end
if LoopCount ~= fix(LoopCount)
    error('������ܵ�����������Ӧ����һ��������');
end
if LoopCount < 1
    error('������ܵ�����������Ӧ����һ�����ڻ���� 1 ��������');
end   

if nargout ~= 2
    error('��������ĸ�������')
end

% 2. ��ʼ��
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

% 3. ������ʼ��Ⱥ
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc); 

[~,row] = min(Solution(1:PopSize,SearchDimension+1));  
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);

% 4. ��ѭ��
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
    % 4.1 ��ʾ��ǰ�����Ӧ�� 
    BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Best fitness:'...
        num2str(BestFitness-AdaptFunc*100,'%e')];
    disp(str);
    nfes = nfes + PopSize;
    
    % 4.2 ���� DE ����
    F = 0.5;
    CR = 0.9;

    % 4.3 ѡ�����
    % ������Ӧ�ȣ�ȡ��ֵ��ʹ��ԭ����С��ֵ���
    fitness = -Solution(1:PopSize, SearchDimension+1);

    % ��׼����Ӧ�ȣ�Z-score��׼����
    fitness_mean = mean(fitness);
    fitness_std = std(fitness);
    fitness_normalized = (fitness - fitness_mean) / fitness_std;

    % �����¶Ȳ���
    T = 1;

    % ʹ��softmax����ѡ�����
    exp_fitness = exp(fitness_normalized / T);
    selection_prob = exp_fitness / sum(exp_fitness);

    % ����ѡ��������ѡ������Ⱥ
    new_population_indices = randsample(PopSize, PopSize, true, selection_prob);
    new_population = Solution(new_population_indices, :);

    % 4.4 �����������
    nrandI = 5;
    r = zeros(PopSize, nrandI);
    for i = 1:PopSize
        original_index = new_population_indices(i);
        available_indices = setdiff(1:PopSize, original_index);
        r(i,:) = available_indices(randperm(PopSize-1, nrandI));
    end

    % 4.5 �������
    V(1:PopSize,1:SearchDimension) = new_population(:,1:SearchDimension) + F * (Solution(r(:,2),1:SearchDimension) - Solution(r(:,3),1:SearchDimension));

    % 4.6 �߽紦��
    V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + new_population(:,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + new_population(:,1:SearchDimension))./2);

    % 4.7 �������
    U = new_population;
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(:,1:SearchDimension) = I .* V(:,1:SearchDimension) + (~I) .* new_population(:,1:SearchDimension);

    % 4.8 �����½�
    U(:,SearchDimension+1) = cec17_func(U(:,1:SearchDimension)',AdaptFunc); 

    % 4.9 ѡ��
    tmp = (U(:,SearchDimension+1) < new_population(:,SearchDimension+1));
    temp = repmat(tmp,1,SearchDimension+1);
    Solution(1:PopSize,:) = temp.*U + (1-temp).*new_population;

    % 4.10 ������ѽ�
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);

    % 4.11 ��¼��ǰ�����Ӧ��ֵ
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
end

% 5. ���ؽ��
Result = Solution(PopSize+1,:);

end



