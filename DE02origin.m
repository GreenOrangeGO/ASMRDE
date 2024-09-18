function [Result,AdaptFuncValue] = DE02(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)

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

% ����ѭ����ʼǰ����һ�����������������
function diversity = calculateDiversity(population)
    [uniquePopulation, ~, ic] = unique(population, 'rows');
    diversity = size(uniquePopulation, 1);
end

% 4. ��ѭ��
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
    % 4.1 ��ʾ��ǰ�����Ӧ�ȺͶ����� (ÿ5000�ε������һ��)
    if mod(k, 5000) == 0 || k == 1 || k == LoopCount
        BestFitness = min(Solution(1:PopSize,SearchDimension+1));
        diversity = calculateDiversity(Solution(1:PopSize,1:SearchDimension));
        str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
            ' Best fitness:' num2str(BestFitness-AdaptFunc*100,'%e')...
            ' Population diversity:' num2str(diversity)];
        disp(str);
    end
    nfes = nfes + PopSize;
    
    % 4.2 ���� DE ����
    F = 0.5;
    CR = 0.9;

    % �����¶Ȳ���
    T = 60.1-0.006*k;

    %500��ѡ����������һ��
    if mod(k,500) ==0 
    
        %��ʼ���浵B
        B = [];

        stagnantCounter = zeros(PopSize, 1);

        % 4.3 ѡ�����
        % ������Ӧ�ȣ�ֱ��ʹ��ԭʼ��Ӧ��ֵ����С��ֵ���ã�
        fitness = -Solution(1:PopSize, SearchDimension+1);

        % ��׼����Ӧ�ȣ�Z-score��׼��������Ҫ�����׼��Ϊ0�����
        fitness_mean = mean(fitness);
        fitness_std = std(fitness);

        if fitness_std == 0
            % �����׼��Ϊ0��˵��������Ӧ��ֵ��ͬ
            % ����������£����ǿ��Ը�ÿ��������ȵ�ѡ�����
            selection_prob = ones(PopSize, 1) / PopSize;
        else
            % �������б�׼���ͺ�������
            fitness_normalized = (fitness - fitness_mean) / fitness_std;
            
            % ʹ��softmax����ѡ�����
            exp_fitness = exp(fitness_normalized/ T);
            selection_prob = exp_fitness / sum(exp_fitness);
            
            % ���С������ �� ��ȷ��û�������
            epsilon = 1e-10;
            selection_prob = selection_prob + epsilon;
            selection_prob = selection_prob / sum(selection_prob);
        end

        % ʹ�� randsample ����
        new_population_indices = randsample(PopSize, PopSize, true, selection_prob);
        new_population = Solution(new_population_indices, :);

        % ����solution��δ��ѡ�и����index
        unselected_indices = setdiff(1:PopSize, new_population_indices);
        
        % ��δ��ѡ�еĸ�����ӵ��浵B��
        B = [B; Solution(unselected_indices, :)];

        % ���B�Ĵ�С����PopSize��������õ�PopSize������
        if size(B, 1) > PopSize
            [~, sorted_indices] = sort(B(:, end));
            B = B(sorted_indices(1:PopSize), :);
        end

        % 4.4 �����������
        nrandI = 5;
        r = zeros(PopSize, nrandI);
        for i = 1:PopSize
            original_index = new_population_indices(i);
            available_indices = setdiff(1:PopSize, original_index);
            r(i,:) = available_indices(randperm(PopSize-1, nrandI));
        end

        % 4.5 �������
        % �Ż���ʹ�þ����������ѭ��
        V = Solution(r(1:PopSize,1),1:SearchDimension) + F * (Solution(r(:,2),1:SearchDimension) - Solution(r(:,3),1:SearchDimension));

        % 4.6 �߽紦��
        % �Ż���ʹ���߼����������������
        lower_bound = repmat(SearchScope(:,1)', [PopSize,1]);
        upper_bound = repmat(SearchScope(:,2)', [PopSize,1]);
        V_pop = new_population(:,1:SearchDimension);

        V_out_lower = V < lower_bound;
        V_out_upper = V > upper_bound;
        V_in_bounds = ~(V_out_lower | V_out_upper);

        V = V .* V_in_bounds + ...
            ((lower_bound + V_pop) / 2) .* V_out_lower + ...
            ((upper_bound + V_pop) / 2) .* V_out_upper;

        % 4.7 �������
        % �Ż���ʹ�þ����������ѭ��
        jRand = repmat(ceil(rand(PopSize,1)*SearchDimension), 1, SearchDimension);
        j = repmat(1:SearchDimension, PopSize, 1);
        I = (rand(PopSize,SearchDimension) < CR) | (j == jRand);
        U(:,1:SearchDimension) = I .* V + (~I) .* new_population(:,1:SearchDimension);

        % 4.8 �����½�
        U(:,SearchDimension+1) = cec17_func(U(:,1:SearchDimension)',AdaptFunc); 

        % 4.9 ѡ��
        tmp = (U(:,SearchDimension+1) < new_population(:,SearchDimension+1));
        temp = repmat(tmp,1,SearchDimension+1);
        Solution(1:PopSize,:) = temp.*U + (1-temp).*new_population; 

    else
        
        if k > 500
        % ���stagnantCounter���Ƿ���Ԫ�ش���32
        replace_indices = find(stagnantCounter > 32);

            if ~isempty(replace_indices)
            % �Ӵ浵B�����ѡ�����
            % ѡ�����
            % ������Ӧ�ȣ�ֱ��ʹ��ԭʼ��Ӧ��ֵ����С��ֵ���ã�
            fitness = -B(:, SearchDimension+1);

            % ��׼����Ӧ�ȣ�Z-score��׼��������Ҫ�����׼��Ϊ0�����
            fitness_mean = mean(fitness);
            fitness_std = std(fitness);

                if fitness_std == 0
                    % �����׼��Ϊ0��˵��������Ӧ��ֵ��ͬ
                    % ����������£����ǿ��Ը�ÿ��������ȵ�ѡ�����
                    selection_prob = ones(height(B), 1) / height(B);
                else
                    % �������б�׼���ͺ�������
                    fitness_normalized = (fitness - fitness_mean) / fitness_std;

                    % ʹ��softmax����ѡ�����
                    exp_fitness = exp(fitness_normalized/ T);
                    selection_prob = exp_fitness / sum(exp_fitness);

                    % ���С������ �� ��ȷ��û�������
                    epsilon = 1e-10;
                    selection_prob = selection_prob + epsilon;
                    selection_prob = selection_prob / sum(selection_prob);
                end

                % ʹ�� randsample ����
                replace_population_indices = randsample(height(B), length(replace_indices), true, selection_prob);
                replacement_individuals = B(replace_population_indices, :);

                %�ǳ����ã�ȫ���˻���ͣ�͸����ƺ�û��һ��ë�ļ�ֵ

                % % �����滻�ĸ�����ӵ��浵B��
                % B = [B; Solution(replace_indices, :)];
                % 
                % % ���B�Ĵ�С����PopSize��������õ�PopSize������
                % if size(B, 1) > PopSize
                %     [~, sorted_indices] = sort(B(:, end));
                %     B = B(sorted_indices(1:PopSize), :);
                % end

                % �滻Solution�ж�Ӧ�ĸ���
                Solution(replace_indices, :) = replacement_individuals;

                % ���ñ��滻�����stagnantCounter
                stagnantCounter(replace_indices) = 0;
            end
        end

        %����ѭ��
        nrandI = 5;
        r = zeros(PopSize,nrandI);
        for i = 1:PopSize
            available_index = [1:i-1 i+1:PopSize];
            sequence_index = randperm (numel(available_index));
            r(i,:) = available_index(sequence_index(1:nrandI));
        end
    
        V(1:PopSize,1:SearchDimension) = Solution(r(1:PopSize,1),1:SearchDimension) + F *(Solution(r(1:PopSize,2),1:SearchDimension) - Solution(r(1:PopSize,3),1:SearchDimension));
    
        % �߽�����
        V(1:PopSize,1:SearchDimension) = ...
         ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
         (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
         (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

        % ����Ϊ���ھ�������Ľ������
        U = Solution(1:PopSize,:);
        jRand = ceil(rand(PopSize,1)*SearchDimension);
        jRand = repmat(jRand,[1,SearchDimension]);
        j = 1:SearchDimension;
        j = repmat(j,[PopSize,1]);
        I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
        U(I) = V(I);
        
        U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc); 
        tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % ���� tmp ��Ϊ�Ӵ��滻�����ı��
        
        if k > 500
            lastGeneration = tmp;
            % Update stagnantCounter based on conditions
            stagnantCounter = stagnantCounter + ((lastGeneration == 0 & tmp == 0) | (lastGeneration == 1 & tmp == 0));
            stagnantCounter(lastGeneration == 0 & tmp == 1) = 0;
        end

        temp = repmat(tmp,1,SearchDimension+1);                                                % tmp ��Dά���ư汾
        Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % �Ӵ��ɹ��ģ�����U���Ӵ����ɹ��ģ���������
    end

    % % 4.1 ��ʾ��ǰ�����Ӧ�ȺͶ����� (ÿ5000�ε������һ��)
    % if mod(k, 5000) == 0 || k == 1 || k == LoopCount
    %     BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    %     diversity = calculateDiversity(Solution(1:PopSize,1:SearchDimension));
    %     str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
    %         ' Best fitness:' num2str(BestFitness-AdaptFunc*100,'%e')...
    %         ' Population diversity:' num2str(diversity)];
    %     disp(str);
    % end
    
    % 4.10 ������ѽ�
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);

    % 4.11 ��¼��ǰ�����Ӧ��ֵ
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
end

% 5. ���ؽ��
Result = Solution(PopSize+1,:);

end