function [vecLdc,vecLdt,wA,ps,pairs,ldcRdm,ldtRdm] = fisherDiscrContrast(xA,yA,xB,yB,nCond)

% [vecLdc,ldcRdm,ldtRdm] = fisherDiscrContrast(xA,yA,xB,yB,nCond)
%
% This function calculates the Fisher Linear Discriminant Contrast and t-test,
% as described in Walther, et al. 2016 (doi:10.1016/j.neuroimage.2015.12.012).
% The code is largely modified from Nili's RSA toolbox, but this version
% outputs both contrasts and t-values, and also doesn't fold over the RDMs,
% which is useful when permuting possible data combinations.
%
% input:
%   xA & xB: design matrices for two sets of data (used to cross-validate)
%   yA & yB: datasets used for cross-validation (must have same nVox)
%   nCond (optional): number of conditions of interest in design matrices
%
% output:
%   vecLdc: vector of pairwise linear discriminant contrast values
%   ldcRdm: lower triangle RDM of contrast values
%   ldtRdm: lower triangle RDM of t-values (st error based on nTimepoints)
%
% dependencies:
%   covdiag: from Nili's RSA toolbox
%            (http://dx.doi.org/10.1371/journal.pcbi.1003553)
%
% ----------------------------
% Andrew Morgan, Muckli Lab, CCNi, University of Glasgow

%% Construct the condition contrast matrix
[Ca,pairs,rdmMask] = constructContrastMatrix(xA,nCond);
Cb = constructContrastMatrix(xB,nCond);

%% fit model to data set A    
betas = ((xA' * xA) \ xA') * yA; % nCond by nVox
fitErrors = yA - (xA * betas);

%% determine Fisher linear discriminant using data set A
%invSa = inv(covdiag(fitErrors)); % nVox by nVox
wA = Ca' * betas * eye(size(betas,2)); % multivariate normalization of set A

%% project set B onto the Fisher discriminant of set A
yB_wA   = yB * wA';                         % project Yb onto wA (nTime x nPairs)
ebb_was = (xB' * xB) \ (xB' * yB_wA);       % nPred x nPairs
eeb_was = yB_wA - (xB * ebb_was);           % nTime x nPairs
nDFb    = size(yB_wA,1) - size(xB,2);       % deg of freedom
esb_was = diag(eeb_was' * eeb_was) / nDFb;  % error variance (nPairs x 1)

ldc     = diag(Cb' * ebb_was);              % nContrasts by 1

se      = sqrt(esb_was .* diag(Cb' * ((xB' * xB) \ Cb)));  % 1 by nContrasts
ldt     = ldc ./ se;                        % nContrasts by 1

ps      = tcdf(-ldt,nDFb);                  % compute ps

%% put values in RDMs
ldcRdm = nan(nCond);
ldtRdm = ldcRdm;
ldcRdm(rdmMask) = ldc;
ldtRdm(rdmMask) = ldt;

ldcRdm(~rdmMask) = 0;
ldtRdm(~rdmMask) = 0;
vecLdt = ldt;
vecLdc = ldc;

