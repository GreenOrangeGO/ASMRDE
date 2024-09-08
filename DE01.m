function [Result,AdaptFuncValue] = DE01(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount)
%% �ô���Ϊ��׼DE��DE��

%(1)�������˵����
% PopSize:          ��ʾ��Ⱥ�Ĵ�С���� 1 ��һ��һ�еġ����ڻ���� 6 ��������

% SearchDimension�� ��ʾ���ӵ�ά������ 1 ��һ��һ�еġ����ڻ���� 1 ������

% SearchScope:      ��ʾ�����������������λ�ø�ά��ȡֵ��Χ���� 1 �� ParticleDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ���
%                   ����ά��Ϊ 3 ������Ⱥ�����ʽΪ [x1Min, x1Max; x2Min, x2Max; x3Min, x3Max]��

% AdaptFunc��       ��ʾ��ѡ�����Ӧ�Ⱥ�����

% LoopCount��       ��ʾ�������ܴ������� 1 ��һ��һ�еġ����ڻ���� 1 �������ò�������ȱʡ����ȱʡֵΪ 1000��
                        
%(2)�������˵����
% Result��          ��ʾ���� LoopCount �ε���֮�����õ������Ž⼰������Ӧ����Ӧ�Ⱥ���ֵ����һ�� 1��(ParticleDimension+1) ����������

% AdaptFuncValue��  ��ʾ���� LoopCount �ε���֮��õ��ĸ��ε������õ�������Ӧ�Ⱥ���ֵ����һ�� 1��LoopCount ����������

% �÷���[Result,AdaptFuncValue] = DE(PopSize,SearchDimension,SearchScope,@AdaptFunc,LoopCount);


%ע�⣺���ȱ�֤���ļ���Matlab������·���У�Ȼ��鿴��ص���ʾ��Ϣ��

%�����ˣ���ǿ��
%����ʱ�䣺2014.3.28
%�ο����ף���

%% ��ջ���
% clc           %������������ǣ���� command window �������
% clear         %������������ǣ���� workspace ��ı���

%% ��ʽ����Ѱ��֮ǰ���������������ȱʡֵ���úͺϷ��Լ�顢������������кϷ��Լ��

%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargin<4 ʱ����������******************************************
%********************************************************************************************************

if nargin < 4
    error('��������ĸ�������')
end

%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargin<4 ʱ����������******************************************
%********************************************************************************************************


%*****�������Ĺ��ܣ�����ĸ������������ PopSize,SearchDimension,SearchScope,AdaptFunc ��ǰ���������ĺϷ���************
%********************************************************************************************************************

%-----�������Ĺ��ܣ�ȷ��������� PopSize �� 1 ��һ��һ�еġ����ڻ���� 6 ������-------------   
[row,colum] = size(PopSize);
if row>1 || colum>1
    error('�������Ⱥ��С����Ӧ����һ�� 1 �� 1 �е�������');
end
if PopSize ~=fix (PopSize) % fix��x)Ϊȡ������
    error('�������Ⱥ��С����Ӧ����һ��������');    
end
if PopSize < 6
    error('�������Ⱥ��С����Ӧ����һ�����ڻ���� 6 ��������');    
end
%-----�������Ĺ��ܣ�ȷ��������� PopSize �� 1 ��һ��һ�еġ����ڻ���� 6 ������-------------
    
    
    
%-----�������Ĺ��ܣ�ȷ��������� SearchDimension �� 1 ��һ��һ�еġ����ڻ���� 1 ������-----  
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
%-----�������Ĺ��ܣ�ȷ��������� SearchDimension �� 1 ��һ��һ�еġ����ڻ���� 1 ������-----
      
    
    
%-----�������Ĺ��ܣ�ȷ��������� SearchScope �� 1 �� SearchDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ���----- 
[row,colum] = size(SearchScope);
if row~=SearchDimension || colum~=2
    error('����ĸ�άȡֵ��Χ����Ӧ����һ�� %g �� 2 �е�ʵ������',SearchDimension);
end
for d = 1:SearchDimension
    if SearchScope(d,2) <= SearchScope(d,1)
        error('����ĸ�άȡֵ��Χ���󣬵� %g ���еĵڶ�������Ӧ�ô��ڵ�һ�����ݡ�',d);
    end
end
%-----�������Ĺ��ܣ�ȷ��������� SearchScope �� 1 �� SearchDimension �����еġ�ÿ���еڶ������ݴ��ڵ�һ�����ݵľ���-----
    
%*****�������Ĺ��ܣ�����ĸ������������ PopSize,SearchDimension,SearchScope,AdaptFunc ��ǰ���������ĺϷ���************
%********************************************************************************************************************      
    
    

%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargin = 4 ʱ������������� LoopCount ��ȱʡֵ**********************
%************************************************************************************************************

if nargin == 4               
    LoopCount = 1000;  
end

%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargin = 4 ʱ������������� LoopCount ��ȱʡֵ**********************
%************************************************************************************************************



%*****�������Ĺ��ܣ��������һ��������� LoopCount �ĺϷ���*************************************************
%**********************************************************************************************************
   
%-------�������Ĺ��ܣ�ȷ��������� LoopCount �� 1 ��һ��һ�еġ����ڻ���� 1 ������----------------------------- 
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
%-------�������Ĺ��ܣ�ȷ��������� LoopCount �� 1 ��һ��һ�еġ����ڻ���� 1 ������------------------------------     
    
%*****�������Ĺ��ܣ��������һ��������� LoopCount �ĺϷ���*************************************************
%**********************************************************************************************************



%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargout ~= 2 ʱ����������*****************************************
%***********************************************************************************************************

if nargout ~= 2
    error('��������ĸ�������')
end

%*****�������Ĺ��ܣ�����������ĸ�����ΧΪ nargout ~= 2 ʱ����������*****************************************
%***********************************************************************************************************


    
%% ��ʼ����Ⱥ
 
%rand('state',0); % ������������ǣ�ÿ�β����������һ����MATLAB��ʾ���޳��á�
                  % ���磺δ��rand('state',0)ʱ������u2 = rand(3,1)���Σ�
                  % ���ֱ����u2 = 0.4057  0.9355  0.9169��u2 = 0.4103  0.8936  0.0579��u2 = 0.3529  0.8132  0.0099�����������ע��ÿ�������������Σ������Ľ�����ܲ�ͬ��
                  % ������rand('state',0)ʱ������u2 = rand(3,1)���Σ�����������ͬ�����ݣ���ÿ��������rand('state',0)��Ȼ�������������Σ������Ľ��������ͬ��
                  % u2 = 0.9501    0.2311    0.6068��u2 =  0.4860    0.8913
                  % 0.7621��u2 = 0.4565    0.0185    0.8214 ���������

rand('state',sum(100*clock)); 
                  
%*****�������Ĺ��ܣ���ʼ�������⼰������Ӧ����Ӧ�Ⱥ���ֵ���Ӷ���ɽ���� Solution �ĳ�ʼ���*****************************
%**********************************************************************************************************************

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ����� Solution �ĽṹΪ������һ�� (PopSize+1)��(SearchDimension+1) �ľ���
% ����ǰ PopSize �У���� i �У�ǰ SearchDimension �ж�Ӧ�ڵ� i ���⣬�� SearchDimension+1 �ж�Ӧ�ڵ� i �������Ӧ�Ⱥ���ֵ��
% ���ڵ� PopSize+1 �У���ǰ SearchDimension �ж�Ӧ��ȫ�����Ž⣬�� SearchDimension+1 �ж�Ӧ��ȫ�����Ž����Ӧ�Ⱥ���ֵ��
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                  

Solution = rand(PopSize+1,SearchDimension+1);                                                  % ���ȣ�������� Solution ȫ����Ϊ [0,1] �����ڵ������
A = generator1(PopSize,SearchDimension);
Solution(1:PopSize,1:SearchDimension) = A;
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
%     Solution(i,SearchDimension+1) = AdaptFunc(Solution(i,1:SearchDimension));             % ��󣬼�����������Ӧ�Ⱥ���ֵ��ѭ�����㷽ʽ��
% end
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc); 
% Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc);  % ��󣬼�����������Ӧ�Ⱥ���ֵ���������㷽ʽ��

% Ѱ����Ӧ�Ⱥ���ֵ��С�Ľ��ھ��� Solution �е�λ��(����)������ɶԽ���� Solution ���һ�еĸ�ֵ
[~,row] = min(Solution(1:PopSize,SearchDimension+1));  
Solution(PopSize+1,1:SearchDimension+1) = Solution(row,1:SearchDimension+1);


%*****�������Ĺ��ܣ���ʼ�������⼰������Ӧ����Ӧ�Ⱥ���ֵ���Ӷ���ɽ���� Solution �ĳ�ʼ���*********************************
%**************************************************************************************************************************


%% ����Ѱ��
%����һ�����Ĺ��ܣ���ʼ�����ڼ�¼ÿ���������������Ӧ�Ⱥ���ֵ������
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
for k = 1:LoopCount
     
    %*****�������Ĺ��ܣ���ʾ�����Ĵ���*********************************
    %******************************************************************
    
%     disp('----------------------------------------------------------')
%     TempStr = sprintf('�� %g �ε���',k);
%     disp(TempStr);
%     disp('----------------------------------------------------------')     
    BestFitness = min(Solution(1:PopSize,SearchDimension+1));
    str=['Function is F' num2str(AdaptFunc) ' Iterations:' num2str(k)...
        ' Best fitness:'...
        num2str(BestFitness-AdaptFunc*100,'%e')];
    disp(str);
    nfes = nfes + PopSize;
    
    %*****�������Ĺ��ܣ������㷨����ز���***********************************************
    %************************************************************************************
    
    
    %*****����������룬���Ը���  ��ͻ�����ӣ� F �ı仯����*****
    
    F = 0.5;
    
    %*****����������룬���Ը���  ��ͻ�����ӣ� F �ı仯����*****
    
    
    %*****����������룬���Ը��� ��������ʣ� CR �ı仯����*****
    
    CR = 0.9;
    
    %*****����������룬���Ը���  ��������ʣ� CR �ı仯����*****
    %*****����������룬���Ը��� �������ѡ����� mutationStrategy ��ȡֵ***********************************
    
%     mutationStrategy = 1;    % mutationStrategy: ��ʾ DE �㷨�ļ��ֳ���ͻ�����
                             % mutationStrategy = 1��   DE/rand/1
                             % mutationStrategy = 2��   DE/best/1
                             % mutationStrategy = 3��   DE/current-to-best/1
    
    %*****����������룬���Ը��� �������ѡ����� mutationStrategy ��ȡֵ***********************************
    
    
    %*****�������Ĺ��ܣ������㷨����ز���***********************************************
    %************************************************************************************
         
    
    
    
    %***** �������Ĺ��ܣ����б������ *************************************************
    %**********************************************************************************
    
%     for i = 1:PopSize
%         
%         %*****�������Ĺ��ܣ�Ϊÿ���� i �� [1 PopSize] �в��� nrandI ��������ȵ������������ i �Բ���� ******************
%         %**************************************************************************************************************
%         
%         nrandI = 5;
%         aa = [1:i-1 i+1:PopSize];
%         bb = randperm (numel(aa));
%         r = aa(bb(1:nrandI));
%         
%         %*****�������Ĺ��ܣ�Ϊÿ���� i �� [1 PopSize] �в��� nrandI ��������ȵ������������ i �Բ���� ******************
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
%                 error('û����ָ���ı�����ԣ��������趨 mutationStrategy ��ֵ');
%         end
%         % �߽�����
%             V(i,1:SearchDimension) = ((V(i,1:SearchDimension)>=SearchScope(:,1)')&(V(i,1:SearchDimension)<=SearchScope(:,2)')).*V(i,1:SearchDimension) + ...
%             (V(i,1:SearchDimension)<SearchScope(:,1)').*(SearchScope(:,1)') + (V(i,1:SearchDimension)>SearchScope(:,2)').*(SearchScope(:,2)');
%     end
    
    %*****�������Ĺ��ܣ�Ϊÿ���� i �� [1 PopSize] �в��� nrandI ��������ȵ������������ i �Բ���� ******************
    %**************************************************************************************************************
    
    nrandI = 5;
    r = zeros(PopSize,nrandI);
    for i = 1:PopSize
        aa = [1:i-1 i+1:PopSize];
        bb = randperm (numel(aa));
        r(i,:) = aa(bb(1:nrandI));
    end
    
    %*****�������Ĺ��ܣ�Ϊÿ���� i �� [1 PopSize] �в��� nrandI ��������ȵ������������ i �Բ���� ******************
    %**************************************************************************************************************
    
    V(1:PopSize,1:SearchDimension) = Solution(r(1:PopSize,1),1:SearchDimension) + F *(Solution(r(1:PopSize,2),1:SearchDimension) - Solution(r(1:PopSize,3),1:SearchDimension));

    % �߽�����
    V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

    


    %***** �������Ĺ��ܣ����б������ *************************************************
    %**********************************************************************************

    
    
    %***** �������Ĺ��ܣ����н������ *************************************************
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
    
    
    % ����Ϊ���ھ�������Ľ������
    U = Solution(1:PopSize,:);
    jRand = ceil(rand(PopSize,1)*SearchDimension);
    jRand = repmat(jRand,[1,SearchDimension]);
    j = 1:SearchDimension;
    j = repmat(j,[PopSize,1]);
    I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
    U(I) = V(I);
    
    
    %***** �������Ĺ��ܣ����н������ *************************************************
    %**********************************************************************************

    
    
    %***** �������Ĺ��ܣ�����ѡ�����**********************************************************
    %**********************************************************************************
    
    
    %*****�������Ĺ��ܣ�����ѭ�����㷽ʽʵ���½���Ӧ�ȼ��㼰�����***********************
%     for i = 1:PopSize
%         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));
%         if U(i,SearchDimension+1) <= Solution(i,SearchDimension+1)
%             Solution(i,1:SearchDimension+1) = U(i,1:SearchDimension+1);
%         end
%     end
    %*****�������Ĺ��ܣ�����ѭ�����㷽ʽʵ���½���Ӧ�ȼ��㼰�����***********************
    
    %*****�������Ĺ��ܣ����ھ������㷽ʽʵ���½���Ӧ�ȼ��㼰�����***********************
%     U(1:PopSize,SearchDimension+1) = AdaptFunc(U(1:PopSize,1:SearchDimension)',AdaptFunc);
%     tmp = (U(1:PopSize,SearchDimension+1)<Solution(1:PopSize,SearchDimension+1));
%     temp = repmat(tmp,1,SearchDimension+1);
%     Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);


%     for i = 1:PopSize
%         U(i,SearchDimension+1) = AdaptFunc(U(i,1:SearchDimension));
%     end
    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc); 
    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % ���� tmp ��Ϊ�Ӵ��滻�����ı��
    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp ��Dά���ư汾
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % �Ӵ��ɹ��ģ�����U���Ӵ����ɹ��ģ���������
    %*****�������Ĺ��ܣ����ھ������㷽ʽʵ���½���Ӧ�ȼ��㼰�����***********************   
    
    
    % Find the fitness-based current best
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    Solution(PopSize+1,:) = Solution(row,:);
    
    %*****�������Ĺ��ܣ�����ѡ�����**********************************************************
    %**********************************************************************************
  
    
    %*****����������ڼ�¼ÿ��Ѱ�ŵ����õ��������Ӧ�Ⱥ���ֵ*******************************
    %**********************************************************************************
    
    AdaptFuncValue(1,k) = Solution(PopSize+1,SearchDimension+1);
    
    %*****����������ڼ�¼ÿ��Ѱ�ŵ����õ��������Ӧ�Ⱥ���ֵ*******************************
    %**********************************************************************************
    
       
end %ȫ��LoopCount�ε���forѭ��������־

%% Ѱ�Ž�������
    
%*****����������ڼ�¼ LoopCount �ε���֮��õ������Ž��***********************************
%**************************************************************************************

Result = Solution(PopSize+1,:);

%*****����������ڼ�¼ LoopCount �ε���֮��õ������Ž��***********************************
%**************************************************************************************

end



