function [ Lh ] = gen_SLIM_Matrix(A,K,type,reg,kadd,gamma)
%  Community detection by SLIM
% - A    sparse adjacency matrix
% - K     number of communities 
% - type  'slim','slimtau','slimk'
% - reg   regularization reg * (average degree)
% - kadd  special for approxiamtion SLIM
% - gamma decay rate
% Email:ningchenying@gmail.com

n = size(A,1);

switch lower(type)
    case 'slim'
        degh_1 = sparse(1./sum(A,2));
        dd = diag(degh_1);
        Tr = dd*A;
        MM = -exp(-gamma)*Tr+speye(n);
        BB = inv(MM);
        Lh = (BB+(BB'))/2;
        Lh = Lh - diag(diag(Lh));
    case 'slimtau'
        tau = (sum(sum(A))/n/n)*reg;
        A= A+tau;
        degh_1 = sparse(1./sum(A,2));
        dd = diag(degh_1);
        Tr = dd*A;
        MM = -exp(-gamma)*Tr+speye(n);
        BB = inv(MM);
        Lh = (BB+(BB'))/2;
        Lh = Lh - diag(diag(Lh));
    case 'slimk'
        gk = kadd;
        degh_1 = sparse(1./sum(A,2));
        dd = diag(degh_1);
        Tr = dd*A;
        W = sparse(n,n);
        CT = Tr;
        for i=1:gk
            W = W+exp(-gamma*i)*CT;
            CT = CT*Tr;
        end
        Lh = (W+(W'))/2;
        Lh = Lh - diag(diag(Lh));
end

end

