function normMuI = compMMR(CM,n)
% normMUI Computes the normalized mutual information between two clusters
%
% CM   confusion matrix

K = size(CM,1);
P = perms(1:K);

l = size(P,1);

Ret = zeros(1,l);

for i=1:K
    A = CM(i,P(:,i))/n;
    %size(A)
    Ret = Ret+A;
end

normMuI = max(Ret);