function [C, pairs, RDMmask] = contructContrastMatrix(X,nCond)

% [C, pairs, RDMmask] = contructContrastMatrix(X,nCond)
%
% Defines a pairwise contrast matrix for regressors of interest in a GLM
% design matrix (X). For instance, A=C'*B wo

%% preparations
nPred = size(X,2);

if ~exist('nCond','var')
    nCond = nPred;
end

% create a square matrix
RDMmask = true(nCond);

% triangularize it (not includeing the diagonal)
RDMmask = tril(RDMmask,-1);

%% construct a contrast matrix C

% index the mask
[pairs(:,2),pairs(:,1)] = find(RDMmask);
nPairs = size(pairs,1);

% Create contrasts
CC = zeros(nCond,nPairs);
ii = sub2ind(size(CC),pairs(:,2),(1:nPairs)');
jj = sub2ind(size(CC),pairs(:,1),(1:nPairs)');
CC(ii) = -1;
CC(jj) = 1;

% put it in C (fits the full design, but only access the conditions)
C = zeros(nPred,nPairs);
C(1:nCond,:) = CC;