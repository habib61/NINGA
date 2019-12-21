function  varargout=BV(varargin)
% Return the score statistic and h2WLS for vectorized image
% heritabilty estimation in terms of WLS estimator
%
% Usage [TS,h2Img,sigmaPImg] = heterV(res,Z,est);
% res      - OLS residual image
% Z        - design matrix for the auxiliary model
% est      - indicator variable that specifies h2 estimation
%
% TS     - Sore statistic image
% h2Img  - heritability image
% h2Img  - Estimated phenotypic varaince image
%__________________________________________________________________________
%
% 
% 
%__________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% May/2015

fields = {'sy','sx','Z'};
for a = 1:nargin,
    eval(sprintf('%s = varargin{a};',fields{a}))
end

%Calculate the residuals
hat               = eye(size(sx,1))-sx*inv(sx'*sx)*sx'; 
res               = hat*sy;
[~,sigmaA,sigmaE] = heterV(res,Z,1);

VC                = [sigmaE';sigmaA'];
%calculate weights
weights           = 1./(Z*VC);

% %calculate \Sigma^-1_{OLS} by 1-step estimator
% F        = res.^2;
% s        = Z\F;
% %Put negative estimates zero
% s(s<0)   = 0;
% weights  = 1./(Z*s);

%calculate the score function
score    =  sx'*(weights.*sy);

%Calculate information matrix
infor    = (sx.*sx)'*weights;
    
%Calculate Beta values

beta     = score./infor;
    
varargout{1} = beta;              % Vectorized beta values  
end