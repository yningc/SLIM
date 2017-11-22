% PolBlogs
load('polblogs.mat')
M = polblogsgiant;
c = polblogsgianty;
K = 2;
% PolBooks
% load('polbooks.mat')
% M = ND;
% c = NT;
% K = 3;

compErr = @(c,e) compMuI(compCM(c,e,K));    % use mutual info as a measure of error/sim.
compErr2 = @(c,e,n) 1-compMMR(compCM(c,e,K),n);

opts = struct('verbose',false,'reg',0.1,'kadd',8,'kmethod',0);

n = size(M,2);

[e,dT] = SLIM(M, K, 'slim', opts);
err = compErr2(c, e,n)

[e1,dT] = SLIM(M, K, 'slimtau', opts);
err1 = compErr2(c, e1,n)

[e2,dT] = SLIM(M, K, 'slimk', opts);
err2 = compErr2(c, e2,n)

opts = struct('verbose',false,'reg',0.1,'kadd',8,'kmethod',1);

[e3,dT] = SLIM(M, K, 'slim', opts);
err3 = compErr2(c, e3,n)