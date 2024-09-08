function [Result,AdaptFuncValue] = DE(PopSize,SearchDimension,SearchScope,AdaptFunc,LoopCount,pp)
rand('state',sum(100*clock)); 
PopSize_Input = SearchDimension;
NFEmax = LoopCount*PopSize_Input;
Solution = rand(PopSize,SearchDimension+1);                                                  % ���ȣ�������� Solution ȫ����Ϊ [0,1] �����ڵ������

% for i = 1:PopSize
%     Solution(i,SearchDimension+1) = AdaptFunc(Solution(i,1:SearchDimension));             % ��󣬼�����������Ӧ�Ⱥ���ֵ��ѭ�����㷽ʽ��
% end
Solution(1:PopSize,SearchDimension+1) = cec17_func(Solution(1:PopSize,1:SearchDimension)',AdaptFunc)'; 

[~,row] = min(Solution(1:PopSize,SearchDimension+1));  
BEST = Solution(row,1:SearchDimension+1);
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
%% main loop
AdaptFuncValue = zeros(1,LoopCount);
nfes = PopSize;
PP = rand(PopSize,SearchDimension);
SS = zeros(PopSize,SearchDimension);
while nfes <=NFEmax
    
    F = 0.5;

    CR = 0.9;
    
    nrandI = 5;
    r = zeros(PopSize,nrandI);
    for i = 1:PopSize
        aa = [1:i-1 i+1:PopSize];
        bb = randperm (numel(aa));
        r(i,:) = aa(bb(1:nrandI));
    end

    V(1:PopSize,1:SearchDimension) = Solution(r(1:PopSize,1),1:SearchDimension) + F *(Solution(r(1:PopSize,2),1:SearchDimension) - Solution(r(1:PopSize,3),1:SearchDimension));

    % �߽�����
    V(1:PopSize,1:SearchDimension) = ...
     ((V(1:PopSize,1:SearchDimension)>=  repmat(SearchScope(:,1)',[PopSize,1])  )&(V(1:PopSize,1:SearchDimension)<=   repmat(SearchScope(:,2)',[PopSize,1]))   ).*V(1:PopSize,1:SearchDimension) + ...
     (V(1:PopSize,1:SearchDimension)<    repmat(SearchScope(:,1)',[PopSize,1])  ).*((   repmat(SearchScope(:,1)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2) + ...
     (V(1:PopSize,1:SearchDimension)>    repmat(SearchScope(:,2)',[PopSize,1])  ).*((   repmat(SearchScope(:,2)',[PopSize,1])  + Solution(1:PopSize,1:SearchDimension))./2);

    % ����Ϊ���ھ�������Ľ������
    
    U = Solution(1:PopSize,1:SearchDimension);
    [U,V,I] = Crossover_strategy_based_on_historical_success_rate(U,V,PP,CR);
%     jRand = ceil(rand(PopSize,1)*SearchDimension);
%     jRand = repmat(jRand,[1,SearchDimension]);
%     j = 1:SearchDimension;
%     j = repmat(j,[PopSize,1]);
%     I = (rand(PopSize,SearchDimension) < repmat(CR,[1,SearchDimension])) | (j == jRand);
%     U(I) = V(I);

    U(1:PopSize,SearchDimension+1) = cec17_func(U(1:PopSize,1:SearchDimension)',AdaptFunc); 
    tmp = (U(1:PopSize,SearchDimension+1) < Solution(1:PopSize,SearchDimension+1));        % ���� tmp ��Ϊ�Ӵ��滻�����ı��
    for i = 1 : PopSize
        for j = 1 : SearchDimension
            SS(i,j) = SS(i,j) + abs(tmp(i)+I(i,j)+1-3*(tmp(i)*I(i,j))); % |a+b+1-3ab| 1 2 3 4; % |3b-a| 0 1 2 3
        end
        if max(SS(i,:)) - min(SS(i,:)) == 0
            PP(i,:) = rand(1,PopSize);
        else
            PP(i,:) = (SS(i,:)-min(SS(i,:)))/(max(SS(i,:)) - min(SS(i,:)));
        end
    end
    temp = repmat(tmp,1,SearchDimension+1);                                                % tmp ��Dά���ư汾
    Solution(1:PopSize,:) = temp.*U(1:PopSize,:) + (1-temp).*Solution(1:PopSize,:);        % �Ӵ��ɹ��ģ�����U���Ӵ����ɹ��ģ���������
    for i = 1:PopSize
        AdaptFuncValue_NFEmax(1,nn) = (AdaptFuncValue_NFEmax(1,nn-1) > Solution(i,SearchDimension+1)).*Solution(i,SearchDimension+1) +...
            (AdaptFuncValue_NFEmax(1,nn-1) <= Solution(i,SearchDimension+1)).*AdaptFuncValue_NFEmax(1,nn-1);
        nn = nn + 1;
    end
    % Find the fitness-based current best
    [~,row] = min(Solution(1:PopSize,SearchDimension+1));
    BEST = Solution(row,:);
    nfes = nfes + PopSize;
end %ȫ��LoopCount�ε���forѭ��������־
kk = 1:LoopCount;
AdaptFuncValue = AdaptFuncValue_NFEmax(PopSize_Input.*kk);
Result = BEST;
str = ['F' num2str(AdaptFunc) ', run ' num2str(pp) ', best_fitness:' num2str(Result(1,end)-100*AdaptFunc,'%e')];
disp(str);
end
function [U,V,I] = Crossover_strategy_based_on_historical_success_rate(U,V,PP,CR)
[a,b] = size(U);
[m,~] = size(CR);
jRand = ceil(rand(a,1)*b);
    jRand = repmat(jRand,[1,b]);
    j = 1:b;
    j = repmat(j,[a,1]);
    if rand <= 0.3
        I = (PP < repmat(CR,[a/m,b])) | (j == jRand);
    else
        I = (rand(a,b) < repmat(CR,[a/m,b])) | (j == jRand);
    end
    U(I) = V(I);
end