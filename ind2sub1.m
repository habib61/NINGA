function varargout = ind2sub1(siz,ndx)
% Works similarly as ind2sub, except that the
% indices are returned all as a single variable,
% hence does not require a priori knowledge
% of the number of dimensions to correctly
% specify nargout.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2012

siz = double(siz);
n = length(siz);

if length(siz)<=n,
  siz = [siz ones(1,n-length(siz))];
else
  siz = [siz(1:n-1) prod(siz(n:end))];
end
k = [1 cumprod(siz(1:end-1))];
indices = zeros(numel(ndx),n);

for i = n:-1:1,
  vi = rem(ndx-1, k(i)) + 1;         
  vj = (ndx - vi)/k(i) + 1; 
  indices(:,i) = vj; 
  ndx = vi;     
end
varargout{1} = indices;