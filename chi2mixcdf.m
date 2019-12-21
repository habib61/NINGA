function P = chi2mixcdf(X,V)
% CDF of a 50:50 mixture of a Chi^2
% distribution and a point mass.
% 
% Usage:
% P = chi2mixcdf(X,V)
% 
% Inputs:
% X: Chi^2 statistic (can be a vector of many values)
% V: Degrees of freedom
% 
% Outputs:
% P: P-value
%
% _____________________________________
% Habib Ganjgahi and Tom Nichols
% 
% Dec/2015

P = zeros(size(X));
P(X == 0) = 0.5;
P(X > 0)  = chi2cdf(X(X > 0), V)/2 + 0.5;