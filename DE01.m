function [Result,AdaptFuncValue] = DE01(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)
%% 该代码为标准DE（DE）

%(1)输入参数说明：
% PopSize:          表示种群的大小，是 1 个一行一列的、大于或等于 6 的整数。

% SearchDimension： 表示粒子的维数，是 1 个一行一列的、大于或等于 1 整数。

% SearchScope:      表示粒子在运算过程中其位置各维的取值范围，是 1 个 ParticleDimension 行两列的、每行中第二个数据大于第一个数据的矩阵。
%                   对于维数为 3 的粒子群，其格式为 [x1Min, x1Max; x2Min, x2Max; x3Min, x3Max]。

% AdaptFunc：       表示所选择的适应度函数。

% LoopCount：       表示迭代的总次数，是 1 个一行一列的、大于或等于 1 整数。该参数可以缺省，其缺省值为 1000。
                        
%(2)输出参数说明：
% Result：          表示经过 LoopCount 次迭代之后所得到的最优解及其所对应的适应度函数值，是一个 1×(ParticleDimension+1) 的行向量。

% AdaptFuncValue：  表示经过 LoopCount 次迭代之后得到的各次迭代所得的最优适应度函数值，是一个 1×LoopCount 的行向量。

% 用法：[Result,AdaptFuncValue] = DE(PopSize,SearchDimension,SearchScope,@AdaptFunc,LoopCount);


%注意：首先保证该文件在Matlab的搜索路径中，然后查看相关的提示信息。

%编制人：杨强大
%编制时间：2014.3.28
%参考文献：无

%% 清空环境
% clc           %该命令的作用是：清除 command window 里的内容
% clear         %该命令的作用是：清除 workspace 里的变量

%% 正式迭代寻优之前，对输入参数进行缺省值设置和合法性检查、对输出参数进行合法性检查

%*****下面代码的功能：当输入参数的个数范围为 nargin<4 时，报错跳出******************************************
%********************************************************************************************************

if nargin < 4
    error('输入参数的个数错误。')
end

%*****上面代码的功能：当输入参数的个数范围为 nargin<4 时，报错跳出******************************************
%********************************************************************************************************


%*****下面代码的功能：检查四个必填输入参数 PopSize,SearchDimension,SearchScope,AdaptFunc 中前三个参数的合法性************
%********************************************************************************************************************

%-----下面代码的功能：确保输入参数 PopSize 是 1 个一行一列的、大于或等于 6 的整数-------------   
[row,colum] = size(PopSize);
if row>1 || colum>1
    error('输入的种群大小错误，应该是一个 1 行 1 列的整数。');
end
if PopSize ~=fix (PopSize) % fix（x)为取整命令
    error('输入的种群大小错误，应该是一个整数。');    
end
if PopSize < 6
    error('输入的种群大小错误，应该是一个大于或等于 6 的整数。');    
end
%-----上面代码的功能：确保输入参数 PopSize 是 1 个一行一列的、大于或等于 6 的整数-------------
    
    
    
%-----下面代码的功能：确保输入参数 SearchDimension 是 1 个一行一列的、大于或等于 1 的整数-----  
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
%-----上面代码的功能：确保输入参数 SearchDimension 是 1 个一行一列的、大于或等于 1 的整数-----
      
    
    
%-----下面代码的功能：确保输入参数 SearchScope 是 1 个 SearchDimension 行两列的、每行中第二个数据大于第一个数据的矩阵----- 
[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('输入的各维取值范围错误，应该是一个 %g 行 2 列的实数矩阵。',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('输入的各维取值范围错误，第 %g 行中的第二个数据应该大于第一个数据。',d);
    end
end
%-----上面代码的功能：确保输入参数 SearchScope 是 1 个 SearchDimension 行两列的、每行中第二个数据大于第一个数据的矩阵-----
    
%*****上面代码的功能：检查四个必填输入参数 PopSize,SearchDimension,SearchScope,AdaptFunc 中前三个参数的合法性************
%********************************************************************************************************************      
    
    

%*****下面代码的功能：当输入参数的个数范围为 nargin = 4 时，设置输入参数 LoopCount 的缺省值**********************
%************************************************************************************************************

if nargin == 4               
    LoopCount = 1000;  
end

%*****上面代码的功能：当输入参数的个数范围为 nargin = 4 时，设置输入参数 LoopCount 的缺省值**********************
%************************************************************************************************************



%*****下面代码的功能：检查其他一个输入参数 LoopCount 的合法性*************************************************
%**********************************************************************************************************
   
%-------下面代码的功能：确保输入参数 LoopCount 是 1 个一行一列的、大于或等于 1 的整数----------------------------- 
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
%-------上面代码的功能：确保输入参数 LoopCount 是 1 个一行一列的、大于或等于 1 的整数------------------------------     
    
%*****上面代码的功能：检查其他一个输入参数 LoopCount 的合法性*************************************************
%**********************************************************************************************************



%*****下面代码的功能：当输出参数的个数范围为 nargout ~= 2 时，报错跳出*****************************************
%***********************************************************************************************************

if nargout ~= 2
    error('输出参数的个数错误。')
end

%*****上面代码的功能：当输出参数的个数范围为 nargout ~= 2 时，报错跳出*****************************************
%***********************************************************************************************************


    
%% 初始化种群
 
%rand('state',0); % 该命令的作用是：每次产生的随机数一样，MATLAB提示不赞成用。
                  % 例如：未用rand('state',0)时，运行u2 = rand(3,1)三次，
                  % 则会分别产生u2 = 0.4057  0.9355  0.9169，u2 = 0.4103  0.8936  0.0579，u2 = 0.3529  0.8132  0.0099三组随机数（注：每次这样运行三次，产生的结果可能不同）
                  % 当用了rand('state',0)时，运行u2 = rand(3,1)三次，则均会产生相同的数据（即每次先运行rand('state',0)，然后这样运行三次，产生的结果均会相同）
                  % u2 = 0.9501    0.2311    0.6068，u2 =  0.4860    0.8913
                  % 0.7621，u2 = 0.4565    0.0185    0.8214 三组随机数

rand('state',sum(100*clock)); 
                  
%*****下面代码的功能：初始化各个解及其所对应的适应度函数值，从而完成解矩阵 Solution 的初始化搭建*****************************
%**********************************************************************************************************************

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 解矩阵 Solution 的结构为：它是一个 (PopSize+1)×(SearchDimension+1) 的矩阵。
% 对于前 PopSize 行，其第 i 行，前 SearchDimension 列对应于第 i 个解，第 SearchDimension+1 列对应于第 i 个解的适应度函数值；
% 对于第 PopSize+1 行，其前 SearchDimension 列对应于全局最优解，第 SearchDimension+1 列对应于全局最优解的适应度函数值。
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                  

Solution = rand(PopSize+1,SearchDimension+1);                                                  % 首先，将解矩阵 Solution 全部设为 [0,1] 区间内的随机数
A = generator1(PopSize,SearchDimension);
Solution(1:PopSize,1:SearchDimension) = A;
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
%     Solution(i,SearchDimension+1) = AdaptFunc(Solution(i,1:SearchDimension));             % 最后，计算各个解的适应度函数值（循环运算方式）
% end
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc); 
% Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc);  % 最后，计算各个解的适应度函数值（矩阵运算方式）

% 寻找适应度函数值最小的解在矩阵 Solution 中的位置(行数)，并完成对解矩阵 Solution 最后一行的赋值
[~,row] = min(Solution(1:PopSize,SearchDimension+1));  
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);


%*****上面代码的功能：初始化各个解及其所对应的适应度函数值，从而完成解矩阵 Solution 的初始化搭建*********************************
%**************************************************************************************************************************


%% 迭代寻优
%以下一句程序的功能：初始化用于记录每步迭代所得最大适应度函数值的向量
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
     
    %*****下面代码的功能：显示迭代的次数*********************************
    %******************************************************************
    
%     disp('----------------------------------------------------------')
%     TempStr = sprintf('第 %g 次迭代',k);
%     disp(TempStr);
%     disp('----------------------------------------------------------')     
    BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Best fitness:'...
        num2str(BestFitness-AdaptFunc*100,'%e')];
    disp(str);
    nfes = nfes + PopSize;
    
    %*****下面代码的功能：设置算法的相关参数***********************************************
    %************************************************************************************
    
    
    %*****更改下面代码，可以更改  （突变因子） F 的变化规律*****
    
    F = 0.5;
    
    %*****更改上面代码，可以更改  （突变因子） F 的变化规律*****
    
    
    %*****更改下面代码，可以更改 （交叉概率） CR 的变化规律*****
    
    CR = 0.9;
    
    %*****更改上面代码，可以更改  （交叉概率） CR 的变化规律*****
    %*****更改下面代码，可以更改 变异策略选择变量 mutationStrategy 的取值***********************************
    
%     mutationStrategy = 1;    % mutationStrategy: 表示 DE 算法的几种常见突变策略
                             % mutationStrategy = 1：   DE/rand/1
                             % mutationStrategy = 2：   DE/best/1
                             % mutationStrategy = 3：   DE/current-to-best/1
    
    %*****更改上面代码，可以更改 变异策略选择变量 mutationStrategy 的取值***********************************
    
    
    %*****上面代码的功能：设置算法的相关参数***********************************************
    %************************************************************************************
         
    
    
    
    %***** 下面代码的功能：进行变异操作 *************************************************
    %**********************************************************************************
    
%     for i = 1:PopSize
%         
%         %*****下面代码的功能：为每个解 i 在 [1 PopSize] 中产生 nrandI 个互不相等的随机数，且与 i 皆不相等 ******************
%         %**************************************************************************************************************
%         
%         nrandI = 5;
%         aa = [1:i-1 i+1:PopSize];
%         bb = randperm (numel(aa));
%         r = aa(bb(1:nrandI));
%         
%         %*****上面代码的功能：为每个解 i 在 [1 PopSize] 中产生 nrandI 个互不相等的随机数，且与 i 皆不相等 ******************
%         %**************************************************************************************************************
%         
%         
%         switch mutationStrategy
%             case 1
%                 V(i,1:SearchDimension) = Solution(r(1),1:SearchDimension) + F*(Solution(r(2),1:SearchDimension)-Solution(r(3),1:SearchDimension));        % mutationStrategy = 1:    DE/rand/1;
%                 
%             case 2
%                 V(i,1:SearchDimension) = Solution(PopSize+1,1:SearchDimension) + F*(Solution(r(1),1:SearchDimension)-Solution(r(2),1:SearchDimension));   % mutationStrategy = 2:    DE/best/1;
%                 
%             case 3
%                 V(i,1:SearchDimension) = Solution(i,1:SearchDimension) + F*(Solution(PopSize+1,1:SearchDimension)-Solution(i,1:SearchDimension)) ...
%                                          + F*(Solution(r(1),1:SearchDimension)-Solution(Solutionr(2),1:SearchDimension));                                         % mutationStrategy = 3:    DE/current-to-best/1;
%                 
%             otherwise
%                 error('没有所指定的变异策略，请重新设定 mutationStrategy 的值');
%         end
%         % 边界修正
%             V(i,1:SearchDimension) = ((V(i,1:SearchDimension)>=SearchScope(:,1)')&(V(i,1:SearchDimension)<=SearchScope(:,2)')).*V(i,1:SearchDimension) + ...
%             (V(i,1:SearchDimension)<SearchScope(:,1)').*(SearchScope(:,1)') + (V(i,1:SearchDimension)>SearchScope(:,2)').*(SearchScope(:,2)');
%     end
    
    %*****下面代码的功能：为每个解 i 在 [1 PopSize] 中产生 nrandI 个互不相等的随机数，且与 i 皆不相等 ******************
    %**************************************************************************************************************
    
    nrandI = 5;
    r = zeros(PopSize,nrandI);
    for i = 1:PopSize
        aa = [1:i-1 i+1:PopSize];
        bb = randperm (numel(aa));
        r(i,:) = aa(bb(1:nrandI));
    end
    
    %*****上面代码的功能：为每个解 i 在 [1 PopSize] 中产生 nrandI 个互不相等的随机数，且与 i 皆不相等 ******************
    %**************************************************************************************************************
    
    V(1:PopSize,1:SearchDimension) = Solution(r(1:PopSize,1),1:SearchDimension) + F *(Solution(r(1:PopSize,2),1:SearchDimension) - Solution(r(1:PopSize,3),1:SearchDimension));

    % 边界修正
    V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

    


    %***** 上面代码的功能：进行变异操作 *************************************************
    %**********************************************************************************

    
    
    %***** 下面代码的功能：进行交叉操作 *************************************************
    %**********************************************************************************
    
%     U = zeros(PopSize,SearchDimension+1);
%     for i = 1:PopSize
%         jRand = ceil(rand*SearchDimension);
%         for j = 1:SearchDimension
%             if rand<CR || j==jRand
%                 U(i,j) = V(i,j);
%             else
%                 U(i,j) = Solution(i,j);
%             end
%         end
%     end
    
    
    % 以下为基于矩阵运算的交叉操作
    U = Solution(1:PopSize,:);
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);
    
    
    %***** 上面代码的功能：进行交叉操作 *************************************************
    %**********************************************************************************

    
    
    %***** 下面代码的功能：进行选择操作**********************************************************
    %**********************************************************************************
    
    
    %*****下面代码的功能：基于循环运算方式实现新解适应度计算及解更新***********************
%     for i = 1:PopSize
%         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));
%         if U(i,SearchDimension+1) <= Solution(i,SearchDimension+1)
%             Solution(i,1:SearchDimension+1) = U(i,1:SearchDimension+1);
%         end
%     end
    %*****上面代码的功能：基于循环运算方式实现新解适应度计算及解更新***********************
    
    %*****下面代码的功能：基于矩阵运算方式实现新解适应度计算及解更新***********************
%     U(1:PopSize,SearchDimension+1) = AdaptFunc(U(1:PopSize,1:SearchDimension)',AdaptFunc);
%     tmp = (U(1:PopSize,SearchDimension+1)<Solution(1:PopSize,SearchDimension+1));
%     temp = repmat(tmp,1,SearchDimension+1);
%     Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);


%     for i = 1:PopSize
%         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));
%     end
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc); 
    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % 产生 tmp 作为子代替换父代的标记
    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp 的D维复制版本
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % 子代成功的，保留U；子代不成功的，保留父代
    %*****上面代码的功能：基于矩阵运算方式实现新解适应度计算及解更新***********************   
    
    
    % Find the fitness-based current best
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);
    
    %*****上面代码的功能：进行选择操作**********************************************************
    %**********************************************************************************
  
    
    %*****下面代码用于记录每步寻优迭代得到的最大适应度函数值*******************************
    %**********************************************************************************
    
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
    
    %*****上面代码用于记录每步寻优迭代得到的最大适应度函数值*******************************
    %**********************************************************************************
    
       
end %全部LoopCount次迭代for循环结束标志

%% 寻优结果的输出
    
%*****下面代码用于记录 LoopCount 次迭代之后得到的最优结果***********************************
%**************************************************************************************

Result = Solution(PopSize+1,:);

%*****上面代码用于记录 LoopCount 次迭代之后得到的最优结果***********************************
%**************************************************************************************

end



