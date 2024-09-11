function [Result,AdaptFuncValue] = MIDE(PopSize,SearchDimension,SScope,FuncNo,MaxGen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数：NP: 种群个体数
%          SearchDimension:测试函数的变量个数
%          SScope: 规模D×2，各个维度的取值范围([Di,min, Di,max])
%          FuncNo: 测试函数序号
%          MaxGen: 最大迭代次数(当每次迭代时评价次数不一时，不适用，改用最大函数评价次数作为终止条件)

% 输出参数：Result: 规模1×D+1，记录每次运行产生的最佳个体，1：D列为最优解的位置，D+1列为对应函数值
%          AdaptFuncValue：规模1×MaxGen，记录每次运行的函数值，对于标准测试集来说，MaxGen取1e4，产生10000个点用于绘制迭代曲线
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 定义随机数种子
rand('state',sum(100*clock));
% 针对CEC测试PopSize_Input等于问题维度
PopSize_Input = SearchDimension;
% 定义最大函数评价次数
NFEmax = MaxGen*PopSize_Input;
% MIDE算法中的记录停滞个体的存档
POPO = [];
% MIDE算法中的连续不更次数阈值T
T = 0.3*SearchDimension; % Set the maximum number of stagnation of each individual
% MIDE算法中被抛弃个体存档的规模
SA = PopSize;
% MIDE算法中的参数F的变化范围
delta = 0.1;
% CR为MIDE算法的交叉率，尊重mDE算法设置
CR = 0.9;
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
% 初始化每个个体的初始停滞次数
POPS(:,SearchDimension+2) = zeros(PopSize,1);
% CCC为记录最佳个体停滞次数，测试时注释。
% CCC = zeros(PopSize,1);
%%%%******************************主要迭代过程*********************************
while nfes < NFEmax

    % 建立种群并集P∪A
    POPM = [POPS;POPO];

    % 获得并集个体数
    PopMSolution = size(POPM,1);

    % 输出命令行
    [BestFitness, ~] = min(POPM(1:end,SearchDimension+1));
    str = ['F' num2str(FuncNo) ' Iter:' num2str(k)...
        ' Eval:' num2str(nfes) ' Best:'...
        num2str(BestFitness-100*FuncNo,'%e')];
    disp(str);

    % 挑选两个随机数并且不等于当前个体
    nrandI = 2;
    rr = zeros(PopSize,nrandI);
    for i = 1:PopSize
        aa = 1:PopMSolution;
        bb = randperm (numel(aa));
        rr(i,:) = aa(bb(1:nrandI));
    end
    % 挑选两个随机数并且不等于当前个体

    % 计算F取值
    MeanF = 0.9 - 0.8 * (nfes-1)/(NFEmax-1);
    F = unifrnd(MeanF - delta,MeanF + delta);
    % 计算F取值

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

    % 利用公式DE/pbest/1 with external archive产生V
    V(1:PopSize,1:SearchDimension) = Xpbest +...
        F*(POPM(rr(1:PopSize,1),1:SearchDimension) - POPM(rr(1:PopSize,2),1:SearchDimension));

    % 将V中超出边界的维度进行处理，处理方式为取对应X和边界范围的平均数
    V(1:PopSize,1:SearchDimension) = ...
        ((V(1:PopSize,1:SearchDimension)>=repmat(SScope(:,1)',[PopSize,1]))&(V(1:PopSize,1:SearchDimension)<=repmat(SScope(:,2)',[PopSize,1]))).*(V(1:PopSize,1:SearchDimension))+...
        (V(1:PopSize,1:SearchDimension)<repmat(SScope(:,1)',[PopSize,1])).*((repmat(SScope(:,1)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2)+...
        (V(1:PopSize,1:SearchDimension)>repmat(SScope(:,2)',[PopSize,1])).*((repmat(SScope(:,2)',[PopSize,1])+POPS(1:PopSize,1:SearchDimension))./2);

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

    % 记录哪些个体被更新
    tmp = (U(1:PopSize,SearchDimension+1) < POPS(1:PopSize,SearchDimension+1));
    temp = repmat(tmp,1,SearchDimension+1);
    
    % 更新连续不更新次数
    POPS(1:end,SearchDimension+2) = (1-tmp).*(POPS(1:end,SearchDimension+2)+1);

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
    POPM = [POPS;POPO];

    % 找出触发迁移机制的个体
    [INDEX,~]= find((POPS(1:PopSize,SearchDimension+2) == T) == 1);

    % [~,indexxx] = min(POPS(1:PopSize,SearchDimension+1));
    % if ismember(indexxx,INDEX)
    %     CCC(indexxx,1) = CCC(indexxx,1) + 1;
    % end

    %***********************迁移机制**************************
    if ~isempty(INDEX)

        % 将触发迁移机制的个体移入存档
        POPO = [POPO;POPS(INDEX,:)];
        POP = numel(INDEX);
        X = zeros(POP,SearchDimension+1);
        [~,indexbest] = sort(POPM(1:end,1+SearchDimension),'ascend');

        % 更新参数p
        p = 1 - 1*(nfes/NFEmax)^1;
        pNP = round(max(1,PopMSolution*p));
        randindex = ceil(rand(POP,1)*pNP);
        randindex = indexbest(randindex);
        Xpbest = POPM(randindex,1:SearchDimension);

        V1 = zeros(POP,SearchDimension);

        nrandI = 2;
        rrr = zeros(POP,nrandI);
        for i = 1:POP
            aa = 1:PopMSolution;
            bb = randperm (numel(aa));
            rrr(i,:) = aa(bb(1:nrandI));
        end

        % 更新参数F
        MeanF = 0.9 - 0.8 * nfes/NFEmax;
        F = unifrnd(MeanF - delta,MeanF + delta);

        % 进行迁移，不交叉，不选择
        V1(1:POP,1:SearchDimension) = Xpbest + ...
            F*(POPM(rrr(1:POP,1),1:SearchDimension) - POPM(rrr(1:POP,2),1:SearchDimension));

        % 将V1中超出边界的维度进行处理，处理方式为取对应X和边界范围的平均数
        V1 = ...
            ((V1>=repmat(SScope(:,1)',[POP,1]))&(V1<=repmat(SScope(:,2)',[POP,1]))).*(V1)+...
            (V1<repmat(SScope(:,1)',[POP,1])).*((repmat(SScope(:,1)',[POP,1])+POPS(INDEX,1:SearchDimension))./2)+...
            (V1>repmat(SScope(:,2)',[POP,1])).*((repmat(SScope(:,2)',[POP,1])+POPS(INDEX,1:SearchDimension))./2);

        % 评价新产生的个体的适应度函数值，并将原来的个体直接替换
        V1(1:POP,SearchDimension+1) = cec17_func(V1(1:POP,1:SearchDimension)',FuncNo)';
        X(1:POP,1:SearchDimension+1) = V1;

        % 记录函数评价值
        for i = 1 : POP
            SearchProcess(1,nn) = ...
                (SearchProcess(1,nn-1) > X(i,SearchDimension+1)).*X(i,SearchDimension+1) +...
                (SearchProcess(1,nn-1) <= X(i,SearchDimension+1)).*SearchProcess(1,nn-1);
            nn = nn + 1;
        end

        % 更新当前评价次数
        nfes = nfes + POP;

        % 重置停滞次数
        X(1:POP,SearchDimension+2) = 0;

        % 更新种群
        POPS(INDEX,:) = X;
    end
    %***********************迁移机制**************************

    % 存档超出NP处理
    if size(POPO,1) > round(SA)
        [~,indexold] = sort(POPO(:,SearchDimension+1),'ascend');
        POPO = POPO(indexold(1:SA),:);
    end
    
    % 更新合并种群
    POPM = [POPS(1:end,:);POPO];

    % 从合并种群中挑出最好的个体并记录
    [~,row] = min(POPM(1:end,SearchDimension+1));
    BEST = POPM(row,1:SearchDimension+1);

    % 更新迭代次数
    k = k + 1;
end
%%%%******************************主要迭代过程*********************************

% 定义绘制迭代曲线代数间隔
kk = 1:MaxGen;

% 输出结果AdaptFuncValue用于绘制曲线
AdaptFuncValue = SearchProcess(PopSize_Input.*kk);

% 输出结果Result用于计算Mean、SD、Best、Worst值，随后用于统计分析
Result = BEST;

end