function [Result,AdaptFuncValue] = SHADE(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)

% [1. 初始化]
rand('state',sum(100*clock)); 

PopSize_Input = SearchDimension;
PopSize = 100;

NFEmax = LoopCount*PopSize_Input;   

A = [];
H = 100;
M = 0.5*ones(2,H);                    % memory cells 在后面的迭代中需要更新
t = 1;                                % 相当于 memory cells 指针
          
Solution = rand(PopSize+1,SearchDimension+1);                                                

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

% [2. 评估初始种群]
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc)';     
                                                                                                           
[~,row] = min(Solution(1:PopSize,SearchDimension+1));
BEST = Solution(row,1:SearchDimension+1);

% [3. 开始主循环]

nn = 1; 
AdaptFuncValue_NFEmax(1,nn) = Solution(1,SearchDimension+1);
nn = nn + 1;

for ii = 2:PopSize
    if AdaptFuncValue_NFEmax(1,nn-1) > Solution(ii,SearchDimension+1)
        AdaptFuncValue_NFEmax(1,nn) = Solution(ii,SearchDimension+1);
        nn = nn + 1;
    else
        AdaptFuncValue_NFEmax(1,nn) = AdaptFuncValue_NFEmax(1,nn-1);
        nn = nn + 1;
    end
end

nfes = PopSize;
k = 1;
while nfes < NFEmax
     
    BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Evaluation times of consumption:' num2str(nfes) ' Best fitness:'...
        num2str(BestFitness-AdaptFunc*100,'%e')];
    disp(str);
    
    % [4. 生成变异向量]
    Index = max(1,ceil(rand(PopSize,1).*H));         % 随机产生每个个体在记忆单元中的索引

    CR(1:PopSize,1) = M(1,Index(1:PopSize,1))' + 0.1*randn(PopSize,1);
    CR = min(1, max(0, CR));                         % truncated to [0 1]
    
    F(1:PopSize,1) = M(2,Index(1:PopSize,1))' + 0.1 * tan(pi * (rand(PopSize,1) - 0.5));
    for i = 1:PopSize
        while F(i,1) <= 0
            F(i,1) = M(2,Index(i,1)) + 0.1 * tan(pi * (rand(1,1) - 0.5)); % if F(i,1)<=0,it is regenerated
        end 
    end
    F = min(1, F);                                                        % if F(i,1)>1,it is truncated to 1
    
    [~, indBest] = sort(Solution(1:PopSize,SearchDimension+1), 'ascend');
    p = 2/PopSize + rand(PopSize,1)*(0.2 - 2/PopSize);
    pNP = max(round(p*PopSize),2);

    randindex = ceil(rand(PopSize,1).*pNP);
    randindex = max(1,randindex);
    Xpbest = Solution(indBest(randindex),:);
    
    archive = [Solution(1:PopSize,:);A];
    

    r0 = 1:PopSize;
    NP0 = length(r0);
    NP1 = NP0;
    [NP2,~] = size(archive);
    
    r = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : 1000
        pos = (r == r0);
        if sum(pos) == 0
            break;
        else                                      % regenerate r if it is equal to r0
            r(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000                               % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    rr = floor(rand(1, NP0) * NP2) + 1;
    for i = 1 : 1000
        pos = ((rr == r) | (rr == r0));
        if sum(pos)==0
            break;
        else                                       % regenerate rr if it is equal to r0 or r
            rr(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000                                % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
    
     
     V(1:PopSize,1:SearchDimension) = Solution(1:PopSize,1:SearchDimension) + repmat(F(1:PopSize,1),[1,SearchDimension]).*(Xpbest(1:PopSize,1:SearchDimension)-Solution(1:PopSize,1:SearchDimension)) ...
            + repmat(F(1:PopSize,1),[1,SearchDimension]).*(Solution(r(1:PopSize),1:SearchDimension)-archive(rr(1:PopSize),1:SearchDimension));
         
     V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);
    
    % [5. 交叉操作]
    U = Solution(1:PopSize,:);
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);
        
    % [6. 评估新解]
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc)';       
                                                                                                   
    nfes = nfes + PopSize;
    
    % [7. 选择]
    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % 产生 tmp 作为子代替换父代的标记
    SF(1:PopSize,:) = tmp.*F(1:PopSize,:);                                                 % 在列向量 SF 中，如果第 i 个解找到了更好的解，则 SF 第 i 个元素保留当前 F，否则设置为 0
    SCR(1:PopSize,:) = tmp.*CR(1:PopSize,:);                                               % 在列向量 SCR 中，如果第 i 个解找到了更好的解，则 SCR 第 i 个元素保留当前 CR，否则设置为 0
    zeyeta(1:PopSize,1) = tmp.*(abs(Solution(1:PopSize,SearchDimension+1) - U(1:PopSize,SearchDimension+1))); % 记录成功子代与父代的差值
    
    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp 的D维复制版本
    FIndex = find(tmp(1:PopSize,end) == 1);              
    FaSolution = Solution(FIndex(1:end),:);                                                % FaSolution 中保存失败的个体
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % 子代成功的，保留U；子代不成功的，保留父代
    for i = 1:PopSize
        AdaptFuncValue_NFEmax(1,nn) = (AdaptFuncValue_NFEmax(1,nn-1) > Solution(i,SearchDimension+1)).*Solution(i,SearchDimension+1) +...
            (AdaptFuncValue_NFEmax(1,nn-1) <= Solution(i,SearchDimension+1)).*AdaptFuncValue_NFEmax(1,nn-1);
        nn = nn + 1;
    end

    % [8. 更新存档和参数]
    A = [A;FaSolution];
    [~, IX]= unique(A, 'rows');
    if length(IX) < size(A, 1)                          % There exist some duplicate solutions
        A = A(IX, :);
    end
    
    if size(A, 1) > PopSize                             % randomly remove some solutions
        rndpos = randperm(size(A, 1)); 
        rndpos = rndpos(1 : PopSize);
        A = A(rndpos, :);
    end
    
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);
   
    Mnew = M;
    if  sum(tmp) == 0                           % 这里用 0，意味着这里的 SF、SCR 是空集，就是没有一个成功的参数       
        Mnew(1,t) = M(1,t);
        Mnew(2,t) = M(2,t); 
    else       
        w(1:PopSize,1) = zeyeta(1:PopSize,1)/sum(zeyeta);                                      % 计算权重 w
        mSCR = sum(w(1:PopSize,1).*SCR(1:PopSize,1));                                          % 求出 mSCR 
        mSF = sum(w(1:PopSize,1).*SF(1:PopSize,1).^2)/sum(w(1:PopSize,1).*SF(1:PopSize,1));    % 求出 mSF       
        Mnew(1,t) = mSCR;
        Mnew(2,t) = mSF;
        t = t + 1;
        if t > H
            t = 1;
        end
    end
    M = Mnew;
    
    % [9. 检查终止条件]
    if nfes <= NFEmax
        [~,row] = min(Solution(1:end,SearchDimension+1));
        BEST = Solution(row,1:SearchDimension+1);
    else
        [~,row] = min(Solution(1:end,SearchDimension+1));
        BEST = Solution(row,1:SearchDimension+1);
    end
end % End while

% [10. 输出结果]
kk = 1:LoopCount;
AdaptFuncValue = AdaptFuncValue_NFEmax(PopSize_Input.*kk);

Result = BEST;

end