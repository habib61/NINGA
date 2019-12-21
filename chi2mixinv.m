function X = chi2mixinv(P,V)
% Inverse CDF of a 50:50 mixture of a Chi^2
% distribution and a point mass.
%
% Usage:
% X = chi2mixcdf(P,V)
%
% Inputs:
% P: P-value (can be a vector of many values)
% V: Degrees of freedom
%
% Outputs:
% X: Chi^2 statistic
%
% _____________________________________
% Habib Ganjgahi and Tom Nichols
% 
% Dec/2015

tonan = P < 0 | P > 1;
if any(tonan),
    warning('Some P-values are not between 0 and 1. The statistic will be marked as NaN.')
end

pm = P <= 0.5;
X = zeros(size(P));
X( pm) = chi2inv(0,V);
X(~pm) = chi2inv(2*(P(~pm)-0.5),V);
X(tonan) = NaN;
