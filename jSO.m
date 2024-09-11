function [Result,AdaptFuncValue] = jSO(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)
%% 该代码为 jSO

% (1)输入参数说明：
% PopSize:          表示种群的大小，是 1 个一行一列的、大于或等于 6 的整数。

% SearchDimension： 表示粒子的维数，是 1 个一行一列的、大于或等于 1 整数。

% SearchScope:      表示粒子在运算过程中其位置各维的取值范围，是 1 个 ParticleDimension 行两列的、每行中第二个数据大于第一个数据的矩阵。
%                   对于维数为 3 的粒子群，其格式为 [x1Min, x1Max; x2Min, x2Max; x3Min, x3Max]。

% AdaptFunc：       表示所选择的适应度函数。

% LoopCount：       表示迭代的总次数，是 1 个一行一列的、大于或等于 1 整数。该参数可以缺省，其缺省值为 1000。

% (2)输出参数说明：
% Result：          表示经过 LoopCount 次迭代之后所得到的最优解及其所对应的适应度函数值，是一个 1×(ParticleDimension+1) 的行向量。

% AdaptFuncValue：  表示经过 LoopCount 次迭代之后得到的各次迭代所得的最优适应度函数值，是一个 1×LoopCount 的行向量。

% 用法：[Result,AdaptFuncValue] = jSO(PopSize,SearchDimension,SearchScope,@AdaptFunc,LoopCount);


% 编 制 人：杨强大，刘鹏
% 编制时间：2021.5.13
% 参考文献：Brest Janez, Mau?ec Mirjam Sepesy, Bo?kovi? Borko.（姓后名前，或者说并未改变期刊所谓姓名的先后次序） Single objective real-parameter optimization: algorithm jSO. Proceedings of the 2017 IEEE Congress on Evolutionary Computation, 2017, 1311–1318.

%% 清空环境
% clc           % 该命令的作用是：清除 command window 里的内容
% clear         % 该命令的作用是：清除 workspace 里的变量

%% 正式迭代寻优之前，对输入参数进行缺省值设置和合法性检查、对输出参数进行合法性检查

% *****下面代码的功能：当输入参数的个数范围为 nargin<4 时，报错跳出******************************************
% ********************************************************************************************************

if nargin < 4
    error('输入参数的个数错误。')
end

% ***** 上面代码的功能：当输入参数的个数范围为 nargin<4 时，报错跳出 ******************************************
% ********************************************************************************************************


% ***** 下面代码的功能：检查四个必填输入参数 PopSize,SearchDimension,SearchScope,AdaptFunc 中前三个参数的合法性 ************
% ********************************************************************************************************************

% ----- 下面代码的功能：确保输入参数 PopSize 是 1 个一行一列的、大于或等于 6 的整数 -------------
[row,colum] = size(PopSize);
if row>1 || colum>1
    error('输入的种群大小错误，应该是一个 1 行 1 列的整数。');
end
if PopSize ~= fix (PopSize) % fix（x)为取整命令
    error('输入的种群大小错误，应该是一个整数。');
end
if PopSize < 6
    error('输入的种群大小错误，应该是一个大于或等于 6 的整数。');
end
% ----- 上面代码的功能：确保输入参数 PopSize 是 1 个一行一列的、大于或等于 6 的整数 -------------



% ----- 下面代码的功能：确保输入参数 SearchDimension 是 1 个一行一列的、大于或等于 1 的整数 -----
[row,colum] = size(SearchDimension);
if row>1 || colum>1
    error('输入的搜索空间维数错误，应该是一个 1 行 1 列的整数。');
end
if SearchDimension ~= fix(SearchDimension) % fix（x)为取整命令
    error('输入的搜索空间维数错误，应该是一个整数。');
end
if SearchDimension < 1
    error('输入的搜索空间维数错误，应该是一个大于或等于 1 的整数。');
end
% ----- 上面代码的功能：确保输入参数 SearchDimension 是 1 个一行一列的、大于或等于 1 的整数 -----



% ----- 下面代码的功能：确保输入参数 SearchScope 是 1 个 SearchDimension 行两列的、每行中第二个数据大于第一个数据的矩阵 -----
[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('输入的各维取值范围错误，应该是一个 %g 行 2 列的实数矩阵。',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('输入的各维取值范围错误，第 %g 行中的第二个数据应该大于第一个数据。',d);
    end
end
% ----- 上面代码的功能：确保输入参数 SearchScope 是 1 个 SearchDimension 行两列的、每行中第二个数据大于第一个数据的矩阵 -----

% ***** 上面代码的功能：检查四个必填输入参数 PopSize,SearchDimension,SearchScope,AdaptFunc 中前三个参数的合法性 **********
% ********************************************************************************************************************



% ***** 下面代码的功能：当输入参数的个数范围为 nargin = 4 时，设置输入参数 LoopCount 的缺省值 ********************
% ************************************************************************************************************

if nargin == 4
    LoopCount = 1000;
end

% ***** 上面代码的功能：当输入参数的个数范围为 nargin = 4 时，设置输入参数 LoopCount 的缺省值 ********************
% ************************************************************************************************************



% ***** 下面代码的功能：检查其他一个输入参数 LoopCount 的合法性 ***********************************************
% **********************************************************************************************************

% ------- 下面代码的功能：确保输入参数 LoopCount 是 1 个一行一列的、大于或等于 1 的整数 -------------------------
[row,colum] = size(LoopCount);
if row>1 || colum>1
    error('输入的总迭代次数错误，应该是一个 1 行 1 列的整数。');
end
if LoopCount ~= fix(LoopCount) % fix（x)为取整命令
    error('输入的总迭代次数错误，应该是一个整数。');
end
if LoopCount < 1
    error('输入的总迭代次数错误，应该是一个大于或等于 1 的整数。');
end
% ------- 上面代码的功能：确保输入参数 LoopCount 是 1 个一行一列的、大于或等于 1 的整数 -------------------------

% ***** 上面代码的功能：检查其他一个输入参数 LoopCount 的合法性 ***********************************************
% **********************************************************************************************************



% ***** 下面代码的功能：当输出参数的个数范围为 nargout ~= 2 时，报错跳出 ***************************************
% ***********************************************************************************************************

if nargout ~= 2
    error('输出参数的个数错误。')
end

% ***** 上面代码的功能：当输出参数的个数范围为 nargout ~= 2 时，报错跳出 ****************************************
% ***********************************************************************************************************


%% 初始化种群

% rand('state',0); % 该命令的作用是：每次产生的随机数一样，MATLAB 提示不赞成用。
% 例如：未用 rand('state',0) 时，运行 u2 = rand(3,1) 三次，
% 则会分别产生 u2 = 0.4057  0.9355  0.9169，u2 = 0.4103  0.8936  0.0579，u2 = 0.3529  0.8132  0.0099 三组随机数（注：每次这样运行三次，产生的结果可能不同）
% 当用了 rand('state',0)时，运行u2 = rand(3,1) 三次，则均会产生相同的数据（即每次先运行 rand('state',0)，然后这样运行三次，产生的结果均会相同）
% u2 = 0.9501    0.2311    0.6068，u2 =  0.4860    0.8913    0.7621，u2 = 0.4565    0.0185    0.8214 三组随机数

rand('state',sum(100*clock));

PopSize_Input = PopSize;

% ***** 下面代码的功能：设置算法相关的固定参数 ********************************************************************************
% **************************************************************************************************************************

p_max = 0.25;
p_min = p_max/2;

A = [];
H = 5;
M = [0.8 0.8 0.8 0.8 0.9;
    0.3 0.3 0.3 0.3 0.9];  % memory cells 在后面的迭代中需要更新
t = 1;                                          % 相当于 memory cells 指针

NFEmax = LoopCount*PopSize_Input;               % 相当于根据算法调用时输入的 LoopCount 和 PopSize 计算总适应度评价次数 NFEmax。这么处理的原因是，jSO 算法采取了
% 变种群策略，为了比较的公平，将该算法的终止条件设置为 NFEmax，且使之跟参与比较的算法适应度评价次数相等。

PopSizemax = ceil(25*log(SearchDimension)*sqrt(SearchDimension));
PopSizemin = 4;
PopSize = PopSizemax;                           % 相当于计算完 NFEmax 之后， 给出初始化时 PopSize 的真实大小，即 PopSize = PopSizemax


% ***** 上面代码的功能：设置算法相关的固定参数 ********************************************************************************
% **************************************************************************************************************************


% ***** 下面代码的功能：初始化各个解及其所对应的适应度函数值，从而完成解矩阵 Solution 的初始化搭建 ***************************
% **********************************************************************************************************************

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 解矩阵 Solution 的结构为：它是一个 (PopSize+1)×(SearchDimension+1) 的矩阵。
% 对于前 PopSize 行，其第 i 行，前 SearchDimension 列对应于第 i 个解，第 SearchDimension+1 列对应于第 i 个解的适应度函数值；
% 对于第 PopSize+1 行，其前 SearchDimension 列对应于全局最优解，第 SearchDimension+1 列对应于全局最优解的适应度函数值。
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solution = rand(PopSize+1,SearchDimension+1);                                                  % 首先，将解矩阵 Solution 全部设为 [0,1] 区间内的随机数

for d = 1:SearchDimension
    if SearchScope(d,1)==-inf && SearchScope(d,2)~=inf
        Solution(:,d) = unifrnd(SearchScope(d,2)-10^10,SearchScope(d,2),PopSize,1);            % 然后，当 SearchScope(d,1)==-inf && SearchScope(d,2)~=inf 时，
        % 在 [SearchScope(d,2)-10^10,SearchScope(d,2)] 范围内，随机产生 PopSize 个解的第 d 维
    end

    if SearchScope(d,1)~=-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(SearchScope(d,1),SearchScope(d,1)+10^10,PopSize,1);            % 然后，当 SearchScope(d,1)~=-inf && SearchScope(d,2)==inf 时，
        % 在 [SearchScope(d,1),SearchScope(d,1)+10^10] 范围内，随机产生 PopSize 个解的第 d 维
    end

    if SearchScope(d,1)==-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(-10^10,10^10,PopSize,1);                                       % 然后，当 SearchScope(d,1)==-inf && SearchScope(d,2)==inf 时，
        % 在 [-10^10,10^10] 范围内，随机产生 PopSize 个解的第 d 维
    end

    if SearchScope(d,1)~=-inf && SearchScope(d,2)~=inf
        Solution(:,d) = Solution(:,d)*(SearchScope(d,2)-SearchScope(d,1))+SearchScope(d,1);    % 然后，当 SearchScope(d,1)~=-inf && SearchScope(d,2)~=inf 时，
        % 限定各个解，使其在指定的范围之内
    end
end

% for i = 1:PopSize
%     Solution(i,SearchDimension+1) = AdaptFunc(Solution(i,1:SearchDimension));                            % 最后，计算各个解的适应度函数值（循环运算方式）
% end

% Solution(1:PopSize,SearchDimension+1) = AdaptFunc(Solution(1:PopSize,1:SearchDimension));                % 最后，计算各个解的适应度函数值（矩阵运算方式）
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc)';     % 最后，计算各个解的适应度函数值（矩阵运算方式）

% 寻找适应度函数值最小的解在矩阵 Solution 中的位置(行数)，并完成对解矩阵 Solution 最后一行的赋值

[~,row] = min(Solution(1:PopSize,SearchDimension+1));
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);

% ***** 上面代码的功能：初始化各个解及其所对应的适应度函数值，从而完成解矩阵 Solution 的初始化搭建 *******************************
% **************************************************************************************************************************


%% 迭代寻优
% 以下一句程序的功能：初始化用于记录每代迭代所得最优适应度函数值的向量。需要指出的，对于采取变种群策略的 jSO 算法，这里的代，并不是对应于该算法的代。
% 该算法的总适应度评价次数是 NFEmax，那么我们一共记录 NFEmax 个全局最优适应度函数值，并将其记录于 AdaptFuncValue_NFEmax 中；AdaptFuncValue 中一共记
% 录 LoopCount 个全局最优适应度函数值，它们分别对应于 AdaptFuncValue_NFEmax 中 第 PopSize_Input、2*PopSize_Input、......LoopCount*PopSize_Input 个数据
% AdaptFuncValue = zeros(1,LoopCount);

% 以下一句程序的功能：初始化用于记录 NFEmax 个全局最优适应度函数值的向量（后来取消了该语句，原因是最后一代迭代完成后可能出现 AdaptFuncValue_NFEmax 中记录个数多于 NFEmax 的情况）
% AdaptFuncValue_NFEmax = zeros(1,NFEmax);

% 以下程序的功能：将初始化所得的 PopSize 个全局最优适应度值记录于 AdaptFuncValue_NFEmax 中
nn = 1; % 设置 AdaptFuncValue_NFEmax 指针变量起始值为 1
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


% 以下程序的功能：初始化适应度评价次数计数器 nfes 和 代数计数器 k

nfes = PopSize;
k = 1;

while nfes < NFEmax

    %*****下面代码的功能：显示迭代的代数和适应度评价次数********************
    %********************************************************************

%     disp('----------------------------------------------------------')
%     TempStr = sprintf('第 %g 次迭代',k);
%     disp(TempStr);
%     TempStr1 = sprintf('第 %g 次适应度评价',nfes);
%     disp(TempStr1);
%     disp('----------------------------------------------------------')
    bestfit = min(Solution(1:end,SearchDimension+1));
    str = ['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Evaluation times of consumption:' num2str(nfes) ' Best fitness:'...
        num2str(bestfit-AdaptFunc*100,'%e')];
    disp(str);
    %*****上面代码的功能：显示迭代的代数和适应度评价次数********************
    %********************************************************************


    %*****下面代码的功能：设置算法相关的变化参数********************************************
    %************************************************************************************

    SA = round(2.6*PopSize);                                      % 计算存档 A 的大小

    Index = max(1,ceil(rand(PopSize,1).*H));                      % 随机产生每个个体在记忆单元中的索引

    %*****下面代码的功能：利用正太分布为每个个体产生一个交叉概率 CR *****

    CR = zeros(PopSize,1);

    for i = 1:PopSize
        if M(1,Index(i,1)) == 1000                                        % 这里用 1000 只是代表这个位置的值被 截止 了
            CR(i,1) = 0;                                                   % 在截止的状态下，CR 设为 0，与L_SHADE一致
        else
            CR(i,1) = M(1,Index(i,1)) + 0.1*randn(1,1);
        end
    end

    CR = min(1, max(0, CR));                                              % truncated to [0 1]

    if nfes < 0.25*NFEmax
        CR = max(CR,0.7);
    elseif nfes >= 0.25*NFEmax && nfes < 0.5*NFEmax
        CR = max(CR,0.6);
    end


    %*****上面代码的功能：利用正太分布为每个个体产生一个交叉概率 CR *****


    %*****下面代码的功能：利用柯西分布为每个个体产生一个缩放因子 F *****

    F = zeros(PopSize,1);
    for i = 1:PopSize
        F(i,1) = M(2,Index(i,1)) + 0.1 * tan(pi * (rand(1,1) - 0.5));
        while F(i,1) <= 0
            F(i,1) = M(2,Index(i,1)) + 0.1 * tan(pi * (rand(1,1) - 0.5)); % if F(i,1)<=0,it is regenerated
        end
    end
    F = min(1, F);                                                        % if F(i,1)>1,it is truncated to 1
    
    if nfes < 0.6*NFEmax
        F = min(0.7,F);
    end
    
    
    if nfes < 0.2*NFEmax
        Fw = F.*0.7;
    elseif nfes < 0.4*NFEmax&& nfes>= 0.2*NFEmax
        Fw = F.*1.2;
    else
        Fw = F.*0.8;
    end


    % ***** 上面代码的功能：利用柯西分布为每个个体产生一个缩放因子 F *****


    % ***** 上面代码的功能：设置算法相关的变化参数 ******************************************
    % ************************************************************************************


    % ***** 下面代码的功能：进行变异操作 *******************************************************************************
    % *****************************************************************************************************************

    % *****下面代码的功能：设置 Xpbest *********************************

    p = (p_max-p_min)*nfes/NFEmax + p_min;
    [~, indBest] = sort(Solution(1:PopSize,SearchDimension+1), 'ascend');
    pNP = max(round(p*PopSize),2);

    randindex = ceil(rand(PopSize,1)*pNP);
    randindex = max(1,randindex);
    Xpbest = Solution(indBest(randindex),:);

    % ***** 上面代码的功能：设置 Xpbest *********************************


    archive = [Solution(1:PopSize,:);A];


    % ***** 下面代码的功能：产生索引 r 和 rr  ********************

    r0 = 1:PopSize;
    NP0 = length(r0);
    NP1 = NP0;
    [NP2,~] = size(archive);

    r = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : 1000
        pos = (r == r0);
        if sum(pos) == 0
            break;
        else                                                              % regenerate r1 if it is equal to r0
            r(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000                                                       % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    rr = floor(rand(1, NP0) * NP2) + 1;
    for i = 1 : 1000
        pos = ((rr == r) | (rr == r0));
        if sum(pos)==0
            break;
        else                                                              % regenerate r2 if it is equal to r0 or r1
            rr(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000                                                       % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end


    % ***** 上面代码的功能：产生索引 r 和 rr ********************


    % *****下面代码的功能：产生完全不同的个体，而不仅仅是索引值不同  ********************

    %     V = zeros(PopSize,SearchDimension);
    %     for i = 1:PopSize

    % *****下面代码的功能：为每个解 i 在 Solution(1:PopSize,:) 中随机选择一个 xr1，且 *****
    % *****与 i 不相等 (注意：DE 算法中所谓的不同，通常是指序号不同，而这里是指值不同) ******
    % **********************************************************************************

    %         k1 = (Solution(1:PopSize,:) - repmat(Solution(i,:),[PopSize,1])) == 0;
    %         k1 = sum(k1,2);

    %         FIND = find(k1 == SearchDimension+1);

    %         nrandI = 1;
    %         aa = 1:PopSize;
    %         aa(:,FIND(1:end,1)) = 0;
    %         aa = find(aa);                                   % 这个命令相对于把 aa 中等于 0 的给删除了
    %         if isempty(aa)
    %             r = i;                                       % 如果 aa 是空集，令 r = i
    %         else
    %             bb = randperm (numel(aa));
    %             r = aa(bb(1:nrandI));
    %         end

    % *****上面代码的功能：为每个解 i 在 Solution(1:PopSize,:) 中随机选择一个 xr1，且******
    % *****与 i 不相等 (注意：DE 算法中所谓的不同，通常是指序号不同，而这里是指值不同)*******
    % **********************************************************************************

    % *****下面代码的功能：为每个解 i 在 由 P∪A 构成的 archive 中随机选择一个 xrr，且*****
    % *****与上面的 xr1 及 i 均不相等 (注意：DE 算法中所谓的不同，通常是指序号不同，而******
    % *****这里是指值不同)***************************************************************
    % **********************************************************************************

    %         k2 = (archive(:,1:SearchDimension+1) - repmat(Solution(i,1:SearchDimension+1),[row,1])) == 0;
    %         k2 = sum(k2,2);
    %         FIND = find(k2 == SearchDimension+1);
    %
    %         k3 = (archive(:,1:SearchDimension+1) - repmat(Solution(r(1),1:SearchDimension+1),[row,1])) == 0;
    %         k3 = sum(k3,2);
    %         FINDD = find(k3 == SearchDimension+1);
    %
    %         nrandII = 1;
    %
    %         cc = 1:row;
    %         cc(:,FIND(1:end,1)) = 0;
    %         cc(:,FINDD(1:end,1)) = 0;
    %         cc = find(cc);                                             % 上面 3 个命令，本质上相当于在 archive 中删除了等于 xr1 或 xi 的索引号
    %
    %         if isempty(cc)
    %             rr = i;                                                % 如果 cc 是空集，令 rr = i
    %         else
    %             dd = randperm (numel(cc));
    %             rr = cc(dd(1:nrandII));
    %         end
    %
    % *****上面代码的功能：为每个解 i 在 由 P∪A 构成的 archive 中随机选择一个 xrr，且*****
    % *****与上面的 xr1 及 i 均不相等 (注意：DE 算法中所谓的不同，通常是指序号不同，而******
    % *****这里是指值不同)***************************************************************
    % **********************************************************************************
    %
    %
    %     end


    % ***** 上面代码的功能：产生完全不同的个体，而不仅仅是索引值不同  *******************


    % ***** 下面代码的功能：实现变异操作 *************************************************************************

    V = zeros(PopSize,SearchDimension);   % 因为每一代种群大小都会变化，都会有差的个体从 Solution 中删掉，但是 V 中差的个体没有删掉，
    % 所以如果继续使用上一代的 V 可能会受影响，实际运算确实产生了影响

    V(1:PopSize,1:SearchDimension) = Solution(1:PopSize,1:SearchDimension) +...
        repmat(Fw(1:PopSize,1),[1,SearchDimension]).*(Xpbest(1:PopSize,1:SearchDimension)-Solution(1:PopSize,1:SearchDimension)) ...
        + repmat(F(1:PopSize,1),[1,SearchDimension]).*(Solution(r(1:PopSize),1:SearchDimension)-archive(rr(1:PopSize),1:SearchDimension));

    % ***** 上面代码的功能：实现变异操作 *************************************************************************


    % 以下是将变异向量限制在指定范围内的代码(注意：这里的限制方法与常规的不一样，具体参见相关 Word 文档)

    V(1:PopSize,1:SearchDimension) = ...
        ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
        (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
        (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

    % ***** 上面代码的功能：进行变异操作 ********************************************************************************
    % *****************************************************************************************************************


    % ***** 下面代码的功能：进行交叉操作 *************************************************
    % **********************************************************************************

    % 以下为基于循环运算的交叉操作
    % U = zeros(PopSize,SearchDimension+1);
    % for i = 1:PopSize
    % jRand = ceil(rand*SearchDimension);
    % for j = 1:SearchDimension
    % if rand < CR(i,1) || j == jRand
    % U(i,j) = V(i,j);
    % else
    % U(i,j) = Solution(i,j);
    % end
    % end
    % end

    % 以下为基于矩阵运算的两种交叉操作
    % 操作方法1：本质上 CR 相关操作是基于矩阵的，而 jrand 相关操作还是基于循环的
    %     U = zeros(PopSize,SearchDimension+1);
    %     K1 = rand(PopSize,SearchDimension) < CR(PopSize,1);
    %     jRand = ceil(rand(PopSize,1)*SearchDimension);
    %     for i = 1:PopSize
    %         K1(i,jRand(i,1)) = 1;
    %     end
    %     U(1:PopSize,1:SearchDimension) = V(1:PopSize,1:SearchDimension).*K1 + Solution(1:PopSize,1:SearchDimension).*(1-K1);

    % 操作方法2：本质上 CR 和 jrand 相关操作均是基于矩阵的
    U = Solution(1:PopSize,:);
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);

    % ***** 上面代码的功能：进行交叉操作 *************************************************
    % **********************************************************************************


    % ***** 下面代码的功能：进行选择操作 *********************************************************
    % ******************************************************************************************

    % ***** 下面代码的功能：计算新解适应度值的两种方式（即循环运算方式和矩阵运算方式） ***********************************************************************

    %     for i = 1:PopSize
    %         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));                          % 计算各个解的适应度函数值（矩阵运算方式）
    %     end

    %     U(1:PopSize,SearchDimension+1) = AdaptFunc(U(1:PopSize,1:SearchDimension));              % 计算各个解的适应度函数值（矩阵运算方式）
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc)';       % 计算各个解的适应度函数值（矩阵运算方式）
    %       -- CEC 2013、CEC 2014 和 CEC 2017系列测试函数
    % ***** 上面代码的功能：计算新解适应度值的两种方式（即循环运算方式和矩阵运算方式） ***********************************************************************


    % ***** 下面代码的功能：判断并记录历史适应度值，并且记录评价次数 ***********************************************************************

    for i = 1:PopSize
        AdaptFuncValue_NFEmax(1,nn) = (AdaptFuncValue_NFEmax(1,nn-1) > Solution(i,SearchDimension+1)).*Solution(i,SearchDimension+1) +...
            (AdaptFuncValue_NFEmax(1,nn-1) <= Solution(i,SearchDimension+1)).*AdaptFuncValue_NFEmax(1,nn-1);
        nn = nn + 1;
    end

    % ***** 上面代码的功能：判断并记录历史适应度值，并且记录评价次数 ***********************************************************************

    % ***** 下面代码的功能：基于矩阵运算方式实现 SF、SCR、权重、种群 和 A 的更新 ***********************************************************************

    SCR = zeros(PopSize,1);
    SF = zeros(PopSize,1);
    zeyeta = zeros(PopSize,1);
    w = zeros(PopSize,1);

    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % 产生 tmp 作为子代替换父代的标记
    SF(1:PopSize,:) = tmp.*F(1:PopSize,:);                                                 % 在列向量 SF 中，如果第 i 个解找到了更好的解，则 SF 第 i 个元素保留当前 F，否则设置为 0
    SCR(1:PopSize,:) = tmp.*CR(1:PopSize,:);                                               % 在列向量 SCR 中，如果第 i 个解找到了更好的解，则 SCR 第 i 个元素保留当前 CR，否则设置为 0
    zeyeta(1:PopSize,1) = tmp.*(abs(Solution(1:PopSize,SearchDimension+1) - U(1:PopSize,SearchDimension+1)));  % 记录成功子代与父代的差值

    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp 的D维复制版本
    FIndex = find(tmp(1:PopSize,end) == 1);
    FaSolution = Solution(FIndex(1:end),:);                                                % FaSolution 中保存失败的个体
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % 子代成功的，保留U；子代不成功的，保留父代


    % 以下进行矩阵 A 的更新
    % 这里需要注意的是，在 JADE 所提文献中，只是相当于说首先把 FaSolution 跟 A
    % 合并，完事从 A 中随机删除多余的 size(A, 1)-Popsize 个失败解。
    % 但通过查阅 JADE 原始程序发现，实际上在上述两步操作之间还加入了一步查重操作，这样可以确保 A 中没有重复元素。
    % 本人在程序中，遵照了原始程序的做法。
    A = [A;FaSolution];
    [~, IX]= unique(A, 'rows');
    if length(IX) < size(A, 1)                                            % There exist some duplicate solutions
        A = A(IX, :);
    end


    if size(A, 1) > SA                                                    % randomly remove some solutions
        rndpos = randperm(size(A, 1));
        rndpos = rndpos(1 : SA);
        A = A(rndpos, :);
    end

    % ***** 上面代码的功能：基于矩阵运算方式实现 SF、SCR、权重、种群 和 A 的更新 ***********************************************************************

    % ***** 上面代码的功能：进行选择操作 *********************************************************
    % ******************************************************************************************


    % 以下程序的功能：更新适应度评价次数计数器 nfes 和 代数计数器 k

    nfes = nfes + PopSize;
    k = k + 1;

    % 以上程序的功能：更新适应度评价次数计数器 nfes 和 代数计数器 k

    % 寻找适应度函数值最优的解在矩阵 Solution 中的位置(行数)，并完成对解矩阵 Solution 最后一行的赋值
    % 补充说明：以下程序目的是使最后一代 Solution 的最后一行记录的是截止到第 NFEmax 次适应度评价所得
    % 到的全局最优解及其适应度函数值


    if nfes <= NFEmax
        [~,row] = min(Solution(1:PopSize,SearchDimension+1));
        Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);
    else
        [~,row] = min(Solution(1:nfes-NFEmax,SearchDimension+1));
        Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);
    end


    % ***** 下面代码的功能：更新 memory cells **********************************************************
    % *************************************************************************************************

    Mnew = M;
    if  sum(tmp) == 0                                          % 这里用 0，意味着这里的 mSF、mSCR 是空集，就是没有一个成功的参数
        Mnew(1,t) = M(1,t);
        Mnew(2,t) = M(2,t);
    else
        w = zeyeta(1:PopSize,1)/sum(zeyeta);                                                      % 计算权重 w
        mSCR = sum(w(1:PopSize,1).*SCR(1:PopSize,1).^2)/sum(w(1:PopSize,1).*SCR(1:PopSize,1));    % 求出 mSCR
        mSF = sum(w(1:PopSize,1).*SF(1:PopSize,1).^2)/sum(w(1:PopSize,1).*SF(1:PopSize,1));       % 求出 mSF

        if M(1,t) == 1000 || (max(SCR) ==0 && sum(tmp) > 1)    % 这里用 0，意味着最大的 SCR 是 0，但并不意味是空集
            Mnew(1,t) = 1000;                                  % 1000 代表原文中的 截止 的意思，这里只是个标志
        else
            Mnew(1,t) = (mSCR + Mnew(1,t))/2;
        end
        Mnew(2,t) = (mSF +  Mnew(2,t))/2;
        if t == H
            Mnew(2,t) = 0.9;
            Mnew(1,t) = 0.9;
        end
        t = t + 1;
        if t > H
            t = 1;
        end
    end
    M = Mnew;

    % ***** 上面代码的功能：更新 memory cells **********************************************************
    % *************************************************************************************************


    % ***** 下面代码的功能：计算 PopSize，对 Solution 进行删除操作； 计算 A 的大小，并对 A 进行删除操作 *********
    % *******************************************************************************************************

    [~, reSolution] = sort(Solution(1:PopSize,SearchDimension+1), 'ascend');

    PopSize = round((PopSizemin-PopSizemax)*nfes/NFEmax + PopSizemax);

    Solution(reSolution(PopSize+1:end),:) = [];

    if size(A,1) > round(2.6*PopSize)
        D = randperm(size(A,1));
        A(D(1:size(A,1)-ceil(2.6*PopSize)),:) = [];
    end

    % ***** 上面代码的功能：计算 PopSize，对 Solution 进行删除操作； 计算 A 的大小，并对 A 进行删除操作 **********
    % *******************************************************************************************************


end % while 循环结束标志




% ***** 下面代码的功能：根据 AdaptFuncValue_NFEmax 得到 AdaptFuncValue ***************
% **********************************************************************************

kk = 1:LoopCount;
AdaptFuncValue = AdaptFuncValue_NFEmax(PopSize_Input.*kk);

% ***** 上面代码的功能：根据 AdaptFuncValue_NFEmax 得到 AdaptFuncValue ***************
% **********************************************************************************



%% 寻优结果的输出

% ***** 下面代码用于记录 NFEmax 次适应度函数评价之后得到的最优结果 *************************
% **************************************************************************************

Result = Solution(PopSize+1,:);

% ***** 上面代码用于记录 NFEmax 次适应度函数评价之后得到的最优结果 *************************
% **************************************************************************************
end