function [Result,AdaptFuncValue] = jSO(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)
%% �ô���Ϊ jSO

% (1)�������˵����
% PopSize:          ��ʾ��Ⱥ�Ĵ�С���� 1 ��һ��һ�еġ����ڻ���� 6 ��������

% SearchDimension�� ��ʾ���ӵ�ά������ 1 ��һ��һ�еġ����ڻ���� 1 ������

% SearchScope:      ��ʾ�����������������λ�ø�ά��ȡֵ��Χ���� 1 �� ParticleDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ���
%                   ����ά��Ϊ 3 ������Ⱥ�����ʽΪ [x1Min, x1Max; x2Min, x2Max; x3Min, x3Max]��

% AdaptFunc��       ��ʾ��ѡ�����Ӧ�Ⱥ�����

% LoopCount��       ��ʾ�������ܴ������� 1 ��һ��һ�еġ����ڻ���� 1 �������ò�������ȱʡ����ȱʡֵΪ 1000��

% (2)�������˵����
% Result��          ��ʾ���� LoopCount �ε���֮�����õ������Ž⼰������Ӧ����Ӧ�Ⱥ���ֵ����һ�� 1��(ParticleDimension+1) ����������

% AdaptFuncValue��  ��ʾ���� LoopCount �ε���֮��õ��ĸ��ε������õ�������Ӧ�Ⱥ���ֵ����һ�� 1��LoopCount ����������

% �÷���[Result,AdaptFuncValue] = jSO(PopSize,SearchDimension,SearchScope,@AdaptFunc,LoopCount);


% �� �� �ˣ���ǿ������
% ����ʱ�䣺2021.5.13
% �ο����ף�Brest Janez, Mau?ec Mirjam Sepesy, Bo?kovi? Borko.���պ���ǰ������˵��δ�ı��ڿ���ν�������Ⱥ���� Single objective real-parameter optimization: algorithm jSO. Proceedings of the 2017 IEEE Congress on Evolutionary Computation, 2017, 1311�C1318.

%% ��ջ���
% clc           % ������������ǣ���� command window �������
% clear         % ������������ǣ���� workspace ��ı���

%% ��ʽ����Ѱ��֮ǰ���������������ȱʡֵ���úͺϷ��Լ�顢������������кϷ��Լ��

% *****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargin<4 ʱ����������******************************************
% ********************************************************************************************************

if nargin < 4
    error('��������ĸ�������')
end

% ***** �������Ĺ��ܣ�����������ĸ�����ΧΪ nargin<4 ʱ���������� ******************************************
% ********************************************************************************************************


% ***** �������Ĺ��ܣ�����ĸ������������ PopSize,SearchDimension,SearchScope,AdaptFunc ��ǰ���������ĺϷ��� ************
% ********************************************************************************************************************

% ----- �������Ĺ��ܣ�ȷ��������� PopSize �� 1 ��һ��һ�еġ����ڻ���� 6 ������ -------------
[row,colum] = size(PopSize);
if row>1 || colum>1
    error('�������Ⱥ��С����Ӧ����һ�� 1 �� 1 �е�������');
end
if PopSize ~= fix (PopSize) % fix��x)Ϊȡ������
    error('�������Ⱥ��С����Ӧ����һ��������');
end
if PopSize < 6
    error('�������Ⱥ��С����Ӧ����һ�����ڻ���� 6 ��������');
end
% ----- �������Ĺ��ܣ�ȷ��������� PopSize �� 1 ��һ��һ�еġ����ڻ���� 6 ������ -------------



% ----- �������Ĺ��ܣ�ȷ��������� SearchDimension �� 1 ��һ��һ�еġ����ڻ���� 1 ������ -----
[row,colum] = size(SearchDimension);
if row>1 || colum>1
    error('����������ռ�ά������Ӧ����һ�� 1 �� 1 �е�������');
end
if SearchDimension ~= fix(SearchDimension) % fix��x)Ϊȡ������
    error('����������ռ�ά������Ӧ����һ��������');
end
if SearchDimension < 1
    error('����������ռ�ά������Ӧ����һ�����ڻ���� 1 ��������');
end
% ----- �������Ĺ��ܣ�ȷ��������� SearchDimension �� 1 ��һ��һ�еġ����ڻ���� 1 ������ -----



% ----- �������Ĺ��ܣ�ȷ��������� SearchScope �� 1 �� SearchDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ��� -----
[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('����ĸ�άȡֵ��Χ����Ӧ����һ�� %g �� 2 �е�ʵ������',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('����ĸ�άȡֵ��Χ���󣬵� %g ���еĵڶ�������Ӧ�ô��ڵ�һ�����ݡ�',d);
    end
end
% ----- �������Ĺ��ܣ�ȷ��������� SearchScope �� 1 �� SearchDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ��� -----

% ***** �������Ĺ��ܣ�����ĸ������������ PopSize,SearchDimension,SearchScope,AdaptFunc ��ǰ���������ĺϷ��� **********
% ********************************************************************************************************************



% ***** �������Ĺ��ܣ�����������ĸ�����ΧΪ nargin = 4 ʱ������������� LoopCount ��ȱʡֵ ********************
% ************************************************************************************************************

if nargin == 4
    LoopCount = 1000;
end

% ***** �������Ĺ��ܣ�����������ĸ�����ΧΪ nargin = 4 ʱ������������� LoopCount ��ȱʡֵ ********************
% ************************************************************************************************************



% ***** �������Ĺ��ܣ��������һ��������� LoopCount �ĺϷ��� ***********************************************
% **********************************************************************************************************

% ------- �������Ĺ��ܣ�ȷ��������� LoopCount �� 1 ��һ��һ�еġ����ڻ���� 1 ������ -------------------------
[row,colum] = size(LoopCount);
if row>1 || colum>1
    error('������ܵ�����������Ӧ����һ�� 1 �� 1 �е�������');
end
if LoopCount ~= fix(LoopCount) % fix��x)Ϊȡ������
    error('������ܵ�����������Ӧ����һ��������');
end
if LoopCount < 1
    error('������ܵ�����������Ӧ����һ�����ڻ���� 1 ��������');
end
% ------- �������Ĺ��ܣ�ȷ��������� LoopCount �� 1 ��һ��һ�еġ����ڻ���� 1 ������ -------------------------

% ***** �������Ĺ��ܣ��������һ��������� LoopCount �ĺϷ��� ***********************************************
% **********************************************************************************************************



% ***** �������Ĺ��ܣ�����������ĸ�����ΧΪ nargout ~= 2 ʱ���������� ***************************************
% ***********************************************************************************************************

if nargout ~= 2
    error('��������ĸ�������')
end

% ***** �������Ĺ��ܣ�����������ĸ�����ΧΪ nargout ~= 2 ʱ���������� ****************************************
% ***********************************************************************************************************


%% ��ʼ����Ⱥ

% rand('state',0); % ������������ǣ�ÿ�β����������һ����MATLAB ��ʾ���޳��á�
% ���磺δ�� rand('state',0) ʱ������ u2 = rand(3,1) ���Σ�
% ���ֱ���� u2 = 0.4057  0.9355  0.9169��u2 = 0.4103  0.8936  0.0579��u2 = 0.3529  0.8132  0.0099 �����������ע��ÿ�������������Σ������Ľ�����ܲ�ͬ��
% ������ rand('state',0)ʱ������u2 = rand(3,1) ���Σ�����������ͬ�����ݣ���ÿ�������� rand('state',0)��Ȼ�������������Σ������Ľ��������ͬ��
% u2 = 0.9501    0.2311    0.6068��u2 =  0.4860    0.8913    0.7621��u2 = 0.4565    0.0185    0.8214 ���������

rand('state',sum(100*clock));

PopSize_Input = PopSize;

% ***** �������Ĺ��ܣ������㷨��صĹ̶����� ********************************************************************************
% **************************************************************************************************************************

p_max = 0.25;
p_min = p_max/2;

A = [];
H = 5;
M = [0.8 0.8 0.8 0.8 0.9;
    0.3 0.3 0.3 0.3 0.9];  % memory cells �ں���ĵ�������Ҫ����
t = 1;                                          % �൱�� memory cells ָ��

NFEmax = LoopCount*PopSize_Input;               % �൱�ڸ����㷨����ʱ����� LoopCount �� PopSize ��������Ӧ�����۴��� NFEmax����ô�����ԭ���ǣ�jSO �㷨��ȡ��
% ����Ⱥ���ԣ�Ϊ�˱ȽϵĹ�ƽ�������㷨����ֹ��������Ϊ NFEmax����ʹ֮������Ƚϵ��㷨��Ӧ�����۴�����ȡ�

PopSizemax = ceil(25*log(SearchDimension)*sqrt(SearchDimension));
PopSizemin = 4;
PopSize = PopSizemax;                           % �൱�ڼ����� NFEmax ֮�� ������ʼ��ʱ PopSize ����ʵ��С���� PopSize = PopSizemax


% ***** �������Ĺ��ܣ������㷨��صĹ̶����� ********************************************************************************
% **************************************************************************************************************************


% ***** �������Ĺ��ܣ���ʼ�������⼰������Ӧ����Ӧ�Ⱥ���ֵ���Ӷ���ɽ���� Solution �ĳ�ʼ��� ***************************
% **********************************************************************************************************************

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ����� Solution �ĽṹΪ������һ�� (PopSize+1)��(SearchDimension+1) �ľ���
% ����ǰ PopSize �У���� i �У�ǰ SearchDimension �ж�Ӧ�ڵ� i ���⣬�� SearchDimension+1 �ж�Ӧ�ڵ� i �������Ӧ�Ⱥ���ֵ��
% ���ڵ� PopSize+1 �У���ǰ SearchDimension �ж�Ӧ��ȫ�����Ž⣬�� SearchDimension+1 �ж�Ӧ��ȫ�����Ž����Ӧ�Ⱥ���ֵ��
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solution = rand(PopSize+1,SearchDimension+1);                                                  % ���ȣ�������� Solution ȫ����Ϊ [0,1] �����ڵ������

for d = 1:SearchDimension
    if SearchScope(d,1)==-inf && SearchScope(d,2)~=inf
        Solution(:,d) = unifrnd(SearchScope(d,2)-10^10,SearchScope(d,2),PopSize,1);            % Ȼ�󣬵� SearchScope(d,1)==-inf && SearchScope(d,2)~=inf ʱ��
        % �� [SearchScope(d,2)-10^10,SearchScope(d,2)] ��Χ�ڣ�������� PopSize ����ĵ� d ά
    end

    if SearchScope(d,1)~=-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(SearchScope(d,1),SearchScope(d,1)+10^10,PopSize,1);            % Ȼ�󣬵� SearchScope(d,1)~=-inf && SearchScope(d,2)==inf ʱ��
        % �� [SearchScope(d,1),SearchScope(d,1)+10^10] ��Χ�ڣ�������� PopSize ����ĵ� d ά
    end

    if SearchScope(d,1)==-inf && SearchScope(d,2)==inf
        Solution(:,d) = unifrnd(-10^10,10^10,PopSize,1);                                       % Ȼ�󣬵� SearchScope(d,1)==-inf && SearchScope(d,2)==inf ʱ��
        % �� [-10^10,10^10] ��Χ�ڣ�������� PopSize ����ĵ� d ά
    end

    if SearchScope(d,1)~=-inf && SearchScope(d,2)~=inf
        Solution(:,d) = Solution(:,d)*(SearchScope(d,2)-SearchScope(d,1))+SearchScope(d,1);    % Ȼ�󣬵� SearchScope(d,1)~=-inf && SearchScope(d,2)~=inf ʱ��
        % �޶������⣬ʹ����ָ���ķ�Χ֮��
    end
end

% for i = 1:PopSize
%     Solution(i,SearchDimension+1) = AdaptFunc(Solution(i,1:SearchDimension));                            % ��󣬼�����������Ӧ�Ⱥ���ֵ��ѭ�����㷽ʽ��
% end

% Solution(1:PopSize,SearchDimension+1) = AdaptFunc(Solution(1:PopSize,1:SearchDimension));                % ��󣬼�����������Ӧ�Ⱥ���ֵ���������㷽ʽ��
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc)';     % ��󣬼�����������Ӧ�Ⱥ���ֵ���������㷽ʽ��

% Ѱ����Ӧ�Ⱥ���ֵ��С�Ľ��ھ��� Solution �е�λ��(����)������ɶԽ���� Solution ���һ�еĸ�ֵ

[~,row] = min(Solution(1:PopSize,SearchDimension+1));
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);

% ***** �������Ĺ��ܣ���ʼ�������⼰������Ӧ����Ӧ�Ⱥ���ֵ���Ӷ���ɽ���� Solution �ĳ�ʼ��� *******************************
% **************************************************************************************************************************


%% ����Ѱ��
% ����һ�����Ĺ��ܣ���ʼ�����ڼ�¼ÿ����������������Ӧ�Ⱥ���ֵ����������Ҫָ���ģ����ڲ�ȡ����Ⱥ���Ե� jSO �㷨������Ĵ��������Ƕ�Ӧ�ڸ��㷨�Ĵ���
% ���㷨������Ӧ�����۴����� NFEmax����ô����һ����¼ NFEmax ��ȫ��������Ӧ�Ⱥ���ֵ���������¼�� AdaptFuncValue_NFEmax �У�AdaptFuncValue ��һ����
% ¼ LoopCount ��ȫ��������Ӧ�Ⱥ���ֵ�����Ƿֱ��Ӧ�� AdaptFuncValue_NFEmax �� �� PopSize_Input��2*PopSize_Input��......LoopCount*PopSize_Input ������
% AdaptFuncValue = zeros(1,LoopCount);

% ����һ�����Ĺ��ܣ���ʼ�����ڼ�¼ NFEmax ��ȫ��������Ӧ�Ⱥ���ֵ������������ȡ���˸���䣬ԭ�������һ��������ɺ���ܳ��� AdaptFuncValue_NFEmax �м�¼�������� NFEmax �������
% AdaptFuncValue_NFEmax = zeros(1,NFEmax);

% ���³���Ĺ��ܣ�����ʼ�����õ� PopSize ��ȫ��������Ӧ��ֵ��¼�� AdaptFuncValue_NFEmax ��
nn = 1; % ���� AdaptFuncValue_NFEmax ָ�������ʼֵΪ 1
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


% ���³���Ĺ��ܣ���ʼ����Ӧ�����۴��������� nfes �� ���������� k

nfes = PopSize;
k = 1;

while nfes < NFEmax

    %*****�������Ĺ��ܣ���ʾ�����Ĵ�������Ӧ�����۴���********************
    %********************************************************************

%     disp('----------------------------------------------------------')
%     TempStr = sprintf('�� %g �ε���',k);
%     disp(TempStr);
%     TempStr1 = sprintf('�� %g ����Ӧ������',nfes);
%     disp(TempStr1);
%     disp('----------------------------------------------------------')
    bestfit = min(Solution(1:end,SearchDimension+1));
    str = ['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Evaluation times of consumption:' num2str(nfes) ' Best fitness:'...
        num2str(bestfit-AdaptFunc*100,'%e')];
    disp(str);
    %*****�������Ĺ��ܣ���ʾ�����Ĵ�������Ӧ�����۴���********************
    %********************************************************************


    %*****�������Ĺ��ܣ������㷨��صı仯����********************************************
    %************************************************************************************

    SA = round(2.6*PopSize);                                      % ����浵 A �Ĵ�С

    Index = max(1,ceil(rand(PopSize,1).*H));                      % �������ÿ�������ڼ��䵥Ԫ�е�����

    %*****�������Ĺ��ܣ�������̫�ֲ�Ϊÿ���������һ��������� CR *****

    CR = zeros(PopSize,1);

    for i = 1:PopSize
        if M(1,Index(i,1)) == 1000                                        % ������ 1000 ֻ�Ǵ������λ�õ�ֵ�� ��ֹ ��
            CR(i,1) = 0;                                                   % �ڽ�ֹ��״̬�£�CR ��Ϊ 0����L_SHADEһ��
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


    %*****�������Ĺ��ܣ�������̫�ֲ�Ϊÿ���������һ��������� CR *****


    %*****�������Ĺ��ܣ����ÿ����ֲ�Ϊÿ���������һ���������� F *****

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


    % ***** �������Ĺ��ܣ����ÿ����ֲ�Ϊÿ���������һ���������� F *****


    % ***** �������Ĺ��ܣ������㷨��صı仯���� ******************************************
    % ************************************************************************************


    % ***** �������Ĺ��ܣ����б������ *******************************************************************************
    % *****************************************************************************************************************

    % *****�������Ĺ��ܣ����� Xpbest *********************************

    p = (p_max-p_min)*nfes/NFEmax + p_min;
    [~, indBest] = sort(Solution(1:PopSize,SearchDimension+1), 'ascend');
    pNP = max(round(p*PopSize),2);

    randindex = ceil(rand(PopSize,1)*pNP);
    randindex = max(1,randindex);
    Xpbest = Solution(indBest(randindex),:);

    % ***** �������Ĺ��ܣ����� Xpbest *********************************


    archive = [Solution(1:PopSize,:);A];


    % ***** �������Ĺ��ܣ��������� r �� rr  ********************

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


    % ***** �������Ĺ��ܣ��������� r �� rr ********************


    % *****�������Ĺ��ܣ�������ȫ��ͬ�ĸ��壬��������������ֵ��ͬ  ********************

    %     V = zeros(PopSize,SearchDimension);
    %     for i = 1:PopSize

    % *****�������Ĺ��ܣ�Ϊÿ���� i �� Solution(1:PopSize,:) �����ѡ��һ�� xr1���� *****
    % *****�� i ����� (ע�⣺DE �㷨����ν�Ĳ�ͬ��ͨ����ָ��Ų�ͬ����������ֵָ��ͬ) ******
    % **********************************************************************************

    %         k1 = (Solution(1:PopSize,:) - repmat(Solution(i,:),[PopSize,1])) == 0;
    %         k1 = sum(k1,2);

    %         FIND = find(k1 == SearchDimension+1);

    %         nrandI = 1;
    %         aa = 1:PopSize;
    %         aa(:,FIND(1:end,1)) = 0;
    %         aa = find(aa);                                   % �����������ڰ� aa �е��� 0 �ĸ�ɾ����
    %         if isempty(aa)
    %             r = i;                                       % ��� aa �ǿռ����� r = i
    %         else
    %             bb = randperm (numel(aa));
    %             r = aa(bb(1:nrandI));
    %         end

    % *****�������Ĺ��ܣ�Ϊÿ���� i �� Solution(1:PopSize,:) �����ѡ��һ�� xr1����******
    % *****�� i ����� (ע�⣺DE �㷨����ν�Ĳ�ͬ��ͨ����ָ��Ų�ͬ����������ֵָ��ͬ)*******
    % **********************************************************************************

    % *****�������Ĺ��ܣ�Ϊÿ���� i �� �� P��A ���ɵ� archive �����ѡ��һ�� xrr����*****
    % *****������� xr1 �� i ������� (ע�⣺DE �㷨����ν�Ĳ�ͬ��ͨ����ָ��Ų�ͬ����******
    % *****������ֵָ��ͬ)***************************************************************
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
    %         cc = find(cc);                                             % ���� 3 ������������൱���� archive ��ɾ���˵��� xr1 �� xi ��������
    %
    %         if isempty(cc)
    %             rr = i;                                                % ��� cc �ǿռ����� rr = i
    %         else
    %             dd = randperm (numel(cc));
    %             rr = cc(dd(1:nrandII));
    %         end
    %
    % *****�������Ĺ��ܣ�Ϊÿ���� i �� �� P��A ���ɵ� archive �����ѡ��һ�� xrr����*****
    % *****������� xr1 �� i ������� (ע�⣺DE �㷨����ν�Ĳ�ͬ��ͨ����ָ��Ų�ͬ����******
    % *****������ֵָ��ͬ)***************************************************************
    % **********************************************************************************
    %
    %
    %     end


    % ***** �������Ĺ��ܣ�������ȫ��ͬ�ĸ��壬��������������ֵ��ͬ  *******************


    % ***** �������Ĺ��ܣ�ʵ�ֱ������ *************************************************************************

    V = zeros(PopSize,SearchDimension);   % ��Ϊÿһ����Ⱥ��С����仯�������в�ĸ���� Solution ��ɾ�������� V �в�ĸ���û��ɾ����
    % �����������ʹ����һ���� V ���ܻ���Ӱ�죬ʵ������ȷʵ������Ӱ��

    V(1:PopSize,1:SearchDimension) = Solution(1:PopSize,1:SearchDimension) +...
        repmat(Fw(1:PopSize,1),[1,SearchDimension]).*(Xpbest(1:PopSize,1:SearchDimension)-Solution(1:PopSize,1:SearchDimension)) ...
        + repmat(F(1:PopSize,1),[1,SearchDimension]).*(Solution(r(1:PopSize),1:SearchDimension)-archive(rr(1:PopSize),1:SearchDimension));

    % ***** �������Ĺ��ܣ�ʵ�ֱ������ *************************************************************************


    % �����ǽ���������������ָ����Χ�ڵĴ���(ע�⣺��������Ʒ����볣��Ĳ�һ��������μ���� Word �ĵ�)

    V(1:PopSize,1:SearchDimension) = ...
        ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
        (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
        (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

    % ***** �������Ĺ��ܣ����б������ ********************************************************************************
    % *****************************************************************************************************************


    % ***** �������Ĺ��ܣ����н������ *************************************************
    % **********************************************************************************

    % ����Ϊ����ѭ������Ľ������
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

    % ����Ϊ���ھ�����������ֽ������
    % ��������1�������� CR ��ز����ǻ��ھ���ģ��� jrand ��ز������ǻ���ѭ����
    %     U = zeros(PopSize,SearchDimension+1);
    %     K1 = rand(PopSize,SearchDimension) < CR(PopSize,1);
    %     jRand = ceil(rand(PopSize,1)*SearchDimension);
    %     for i = 1:PopSize
    %         K1(i,jRand(i,1)) = 1;
    %     end
    %     U(1:PopSize,1:SearchDimension) = V(1:PopSize,1:SearchDimension).*K1 + Solution(1:PopSize,1:SearchDimension).*(1-K1);

    % ��������2�������� CR �� jrand ��ز������ǻ��ھ����
    U = Solution(1:PopSize,:);
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);

    % ***** �������Ĺ��ܣ����н������ *************************************************
    % **********************************************************************************


    % ***** �������Ĺ��ܣ�����ѡ����� *********************************************************
    % ******************************************************************************************

    % ***** �������Ĺ��ܣ������½���Ӧ��ֵ�����ַ�ʽ����ѭ�����㷽ʽ�;������㷽ʽ�� ***********************************************************************

    %     for i = 1:PopSize
    %         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));                          % ������������Ӧ�Ⱥ���ֵ���������㷽ʽ��
    %     end

    %     U(1:PopSize,SearchDimension+1) = AdaptFunc(U(1:PopSize,1:SearchDimension));              % ������������Ӧ�Ⱥ���ֵ���������㷽ʽ��
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc)';       % ������������Ӧ�Ⱥ���ֵ���������㷽ʽ��
    %       -- CEC 2013��CEC 2014 �� CEC 2017ϵ�в��Ժ���
    % ***** �������Ĺ��ܣ������½���Ӧ��ֵ�����ַ�ʽ����ѭ�����㷽ʽ�;������㷽ʽ�� ***********************************************************************


    % ***** �������Ĺ��ܣ��жϲ���¼��ʷ��Ӧ��ֵ�����Ҽ�¼���۴��� ***********************************************************************

    for i = 1:PopSize
        AdaptFuncValue_NFEmax(1,nn) = (AdaptFuncValue_NFEmax(1,nn-1) > Solution(i,SearchDimension+1)).*Solution(i,SearchDimension+1) +...
            (AdaptFuncValue_NFEmax(1,nn-1) <= Solution(i,SearchDimension+1)).*AdaptFuncValue_NFEmax(1,nn-1);
        nn = nn + 1;
    end

    % ***** �������Ĺ��ܣ��жϲ���¼��ʷ��Ӧ��ֵ�����Ҽ�¼���۴��� ***********************************************************************

    % ***** �������Ĺ��ܣ����ھ������㷽ʽʵ�� SF��SCR��Ȩ�ء���Ⱥ �� A �ĸ��� ***********************************************************************

    SCR = zeros(PopSize,1);
    SF = zeros(PopSize,1);
    zeyeta = zeros(PopSize,1);
    w = zeros(PopSize,1);

    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % ���� tmp ��Ϊ�Ӵ��滻�����ı��
    SF(1:PopSize,:) = tmp.*F(1:PopSize,:);                                                 % �������� SF �У������ i �����ҵ��˸��õĽ⣬�� SF �� i ��Ԫ�ر�����ǰ F����������Ϊ 0
    SCR(1:PopSize,:) = tmp.*CR(1:PopSize,:);                                               % �������� SCR �У������ i �����ҵ��˸��õĽ⣬�� SCR �� i ��Ԫ�ر�����ǰ CR����������Ϊ 0
    zeyeta(1:PopSize,1) = tmp.*(abs(Solution(1:PopSize,SearchDimension+1) - U(1:PopSize,SearchDimension+1)));  % ��¼�ɹ��Ӵ��븸���Ĳ�ֵ

    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp ��Dά���ư汾
    FIndex = find(tmp(1:PopSize,end) == 1);
    FaSolution = Solution(FIndex(1:end),:);                                                % FaSolution �б���ʧ�ܵĸ���
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % �Ӵ��ɹ��ģ�����U���Ӵ����ɹ��ģ���������


    % ���½��о��� A �ĸ���
    % ������Ҫע����ǣ��� JADE ���������У�ֻ���൱��˵���Ȱ� FaSolution �� A
    % �ϲ������´� A �����ɾ������� size(A, 1)-Popsize ��ʧ�ܽ⡣
    % ��ͨ������ JADE ԭʼ�����֣�ʵ������������������֮�仹������һ�����ز�������������ȷ�� A ��û���ظ�Ԫ�ء�
    % �����ڳ����У�������ԭʼ�����������
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

    % ***** �������Ĺ��ܣ����ھ������㷽ʽʵ�� SF��SCR��Ȩ�ء���Ⱥ �� A �ĸ��� ***********************************************************************

    % ***** �������Ĺ��ܣ�����ѡ����� *********************************************************
    % ******************************************************************************************


    % ���³���Ĺ��ܣ�������Ӧ�����۴��������� nfes �� ���������� k

    nfes = nfes + PopSize;
    k = k + 1;

    % ���ϳ���Ĺ��ܣ�������Ӧ�����۴��������� nfes �� ���������� k

    % Ѱ����Ӧ�Ⱥ���ֵ���ŵĽ��ھ��� Solution �е�λ��(����)������ɶԽ���� Solution ���һ�еĸ�ֵ
    % ����˵�������³���Ŀ����ʹ���һ�� Solution �����һ�м�¼���ǽ�ֹ���� NFEmax ����Ӧ����������
    % ����ȫ�����Ž⼰����Ӧ�Ⱥ���ֵ


    if nfes <= NFEmax
        [~,row] = min(Solution(1:PopSize,SearchDimension+1));
        Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);
    else
        [~,row] = min(Solution(1:nfes-NFEmax,SearchDimension+1));
        Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);
    end


    % ***** �������Ĺ��ܣ����� memory cells **********************************************************
    % *************************************************************************************************

    Mnew = M;
    if  sum(tmp) == 0                                          % ������ 0����ζ������� mSF��mSCR �ǿռ�������û��һ���ɹ��Ĳ���
        Mnew(1,t) = M(1,t);
        Mnew(2,t) = M(2,t);
    else
        w = zeyeta(1:PopSize,1)/sum(zeyeta);                                                      % ����Ȩ�� w
        mSCR = sum(w(1:PopSize,1).*SCR(1:PopSize,1).^2)/sum(w(1:PopSize,1).*SCR(1:PopSize,1));    % ��� mSCR
        mSF = sum(w(1:PopSize,1).*SF(1:PopSize,1).^2)/sum(w(1:PopSize,1).*SF(1:PopSize,1));       % ��� mSF

        if M(1,t) == 1000 || (max(SCR) ==0 && sum(tmp) > 1)    % ������ 0����ζ������ SCR �� 0����������ζ�ǿռ�
            Mnew(1,t) = 1000;                                  % 1000 ����ԭ���е� ��ֹ ����˼������ֻ�Ǹ���־
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

    % ***** �������Ĺ��ܣ����� memory cells **********************************************************
    % *************************************************************************************************


    % ***** �������Ĺ��ܣ����� PopSize���� Solution ����ɾ�������� ���� A �Ĵ�С������ A ����ɾ������ *********
    % *******************************************************************************************************

    [~, reSolution] = sort(Solution(1:PopSize,SearchDimension+1), 'ascend');

    PopSize = round((PopSizemin-PopSizemax)*nfes/NFEmax + PopSizemax);

    Solution(reSolution(PopSize+1:end),:) = [];

    if size(A,1) > round(2.6*PopSize)
        D = randperm(size(A,1));
        A(D(1:size(A,1)-ceil(2.6*PopSize)),:) = [];
    end

    % ***** �������Ĺ��ܣ����� PopSize���� Solution ����ɾ�������� ���� A �Ĵ�С������ A ����ɾ������ **********
    % *******************************************************************************************************


end % while ѭ��������־




% ***** �������Ĺ��ܣ����� AdaptFuncValue_NFEmax �õ� AdaptFuncValue ***************
% **********************************************************************************

kk = 1:LoopCount;
AdaptFuncValue = AdaptFuncValue_NFEmax(PopSize_Input.*kk);

% ***** �������Ĺ��ܣ����� AdaptFuncValue_NFEmax �õ� AdaptFuncValue ***************
% **********************************************************************************



%% Ѱ�Ž�������

% ***** ����������ڼ�¼ NFEmax ����Ӧ�Ⱥ�������֮��õ������Ž�� *************************
% **************************************************************************************

Result = Solution(PopSize+1,:);

% ***** ����������ڼ�¼ NFEmax ����Ӧ�Ⱥ�������֮��õ������Ž�� *************************
% **************************************************************************************
end