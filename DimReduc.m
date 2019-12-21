function varargout=DimReduc(varargin)
% 

%__________________________________________________________________________
%

% _________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% July/2015


fields = {'M','S'};
for a = 1:nargin
   eval(sprintf('%s = varargin{a};',fields{a}))
end

nS = size(M,1);

%EVD
%hat = eye(size(M,1))-M*inv(M'*M)*M';
hat = eye(size(M,1))-M*pinv(M);
hat = (hat+hat')/2;
if ~isempty(S)
    if iscell(S)
        tt  = hat*S{1}*hat;
    else
        tt  = hat*S*hat; 
    end
else
    tt  = hat;
end
tt  = (tt+tt')/2;
%EVD
[vector,value] = eig(tt);

%Project to lower dimension
nC    = max(find(diag(value)<1e-6))+1
vec   = vector(:,nC:nS);
value = value(nC:nS,nC:nS);

varargout{1} = vec;
varargout{2} = value;

end


