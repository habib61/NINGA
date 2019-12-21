function [X,Z,eCm,eCx] = NINGA_partition(M,C)
% Partition a design matrix into regressors of interest and
% nuisance according to a given contrast.
% 
% Usage
% [X,Z] = partition(M,C,meth,Y)
% 
% Inputs:
% M    : Design matrix, to be partitioned.
% C    : Contrast that will define the partitioning.
% 
% Outputs:
% X    : Matrix with regressors of interest.
% Z    : Matrix with regressors of no interest.
% eCm  : Effective contrast, equivalent to the original,
%        for the partitioned model [X Z], and considering
%        all regressors.
% eCx  : Same as above, but considering only X.
%
% References:
% * Guttman I. Linear Models: An Introduction. Wiley,
%   New York, 1982
% * Smith SM, Jenkinson M, Beckmann C, Miller K,
%   Woolrich M. Meaningful design and contrast estimability
%   in FMRI. Neuroimage 2007;34(1):127-36.
% * Ridgway GR. Statistical analysis for longitudinal MR

%   imaging of dementia. PhD thesis. 2009.
% * Winkler AM, Ridgway GR, Webster MG, Smith SM,
%   Nichols TE. Permutation inference for the general
%   linear model. Neuroimage 2014 (in press).
% _____________________________________


C2  = null(C');
Q   = pinv(M'*M);
F1  = pinv(C'*Q*C);
Pc  = C*pinv(C'*Q*C)*C'*Q;
C3  = C2 - Pc*C2;
F3  = pinv(C3'*Q*C3);
X   = M*Q*C*F1;
Z   = M*Q*C3*F3;
eCm = vertcat(eye(size(X,2)),...
    zeros(size(Z,2),size(X,2)));
eCx = eCm(1:size(X,2),:);
end