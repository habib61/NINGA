function varargout=VarB(varargin)
% This function returns variance of fixed effect under either mixed effect or fixed effect model.

%__________________________________________________________________________
%

% _________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% December/2014


fields = {'sy','sx','Z'};
for a = 1:nargin,
   eval(sprintf('%s = varargin{a};',fields{a}))
end
%Calculate the full model residuals
hat                 = eye(size(sy,1))-sx*inv(sx'*sx)*sx';
res                 = hat*sy;
nS                  = size(res,1);
if ~isempty(Z)
    %Calculate variance of parameter estimate under the full model
    [~,sigmaA,sigmaE] = heterV(res,Z,1);
    VC                = [sigmaE';sigmaA'];
    %calculate Sigma^-1
    weights           = 1./(Z*VC);
    %Calculate information matrix
    infor             = (sx.*sx)'*weights;
    Sigma             =  1./infor;
else
    Sigma             = mean(res.*res)./(sx'*sx);
end
varargout{1} = Sigma;