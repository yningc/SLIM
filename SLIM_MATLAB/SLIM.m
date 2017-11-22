function [ e, T ] = SLIM(A,K,type,varargin)
%  Community detection by SLIM
% - A    sparse adjacency matrix
% - K     number of communities 
% - type  'ilm' Inverse Laplacian Matrix Method
% Email:ningchenying@gmail.com

options = struct('verbose',false,'perturb',false,...
                'reg',0, ...
                'kadd',0, ...
                'gamma',0.25,...
                'kmethod',0,...
                'eigmethod',0,...
                'normalize',false,'itrNum',5, 'degPert', 0.01);  %default options
if nargin > 3
    % process options
    optNames = fieldnames(options);
    passedOpt = varargin{end};
  
    if ~isstruct(passedOpt), 
        error('Last argument should be an options struct.')
    end
    for fi = reshape(fieldnames(passedOpt),1,[])
        if any(strmatch(fi,optNames))
            options.(fi{:}) = passedOpt.(fi{:});
        else
            error('Unrecognized option name: %s', fi{:})
        end
    end
end

n = size(A,1);

tic

Lh = gen_SLIM_Matrix(A,K,type,options.reg,options.kadd,options.gamma);

fun = @(x) Lh*x;
        
opts.issym = 1;
if options.verbose
    opts.disp = 2;
end

% spectral methods
if options.eigmethod == 0
    [U, ~] = eigs(fun,n,K,'LM',opts);
else
    [U, ~] = eig(full(Lh));
    U = U(:,n:-1:n-K+1);
end

if options.verbose
    kmopts = statset('Display','iter'); 
else
    kmopts = statset('Display','off'); 
end
% clustering methods
if options.kmethod==0 % kmeans
    kmIDX = kmeans(U(:,2:K),K,'replicates',10,...
    'onlinephase','off','Options',kmopts);
    e = reshape(repmat(kmIDX,1,1)',n,[]);
elseif options.kmethod==1 % kmedoids
    kmIDX = kmedoids(U(:,2:K),K,'Algorithm','pam','Distance','euclidean','replicates',10,...
    'onlinephase','off','Options',kmopts);
    e = reshape(repmat(kmIDX,1,1)',n,[]);
elseif options.kmethod==2 % kmeans2
    e = mykmeans(U(:,2:K),K);
end

T = toc;
end

