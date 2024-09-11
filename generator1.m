function A = generator1(NP,D)
if nargin ~=2
    NP = 50;
    D = 50;
end
n = D+1;
j = 2;
q = n-1;
A = oa_permut(q,n,j);

[m,n] = size(A);
if size(A,1)>NP
    A(randperm(m,m-NP),:) = [];
end
if size(A,2)>D
    A(:,randperm(n,n-D)) = [];
end
MAX = max(max(A));
% MIN = min(min(A));
A = A/MAX;
end