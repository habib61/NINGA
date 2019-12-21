function  varargout=heterV(varargin)
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
% December/2014


fields = {'res','Z','est'};
for a = 1:nargin,
    eval(sprintf('%s = varargin{a};',fields{a}))
end
nS      = size(res,1);
Ff      = res.^2;
SigmaP  = mean(Ff);
F       = (Ff./repmat(SigmaP,nS,1))-1;

b       = inv(Z'*Z);
a       = F'*Z;
%c=f'z*inv(z'z)
c       = a*b;
tt      = c(:,1).*a(:,1);
tt2     = c(:,2).*a(:,2);
%calculate T_S statistic image
Ts      = (tt+tt2)/2;
Score   = .5*((Z(:,2)'*F)./SigmaP);
Ts(Score<0)=0;
% %%%estimating H2WLS
if exist('est','var')
    s       =Z\Ff;
    %new modification based on TAOS
    %s(Score<0)=0;
    s(s<0)  = 0;
    weights = 1./((Z*s).^2);
    a       = sum(weights);
    b       = Z(:,2)'*weights;
    c       = Z(:,2)'.^2*weights;
    d       = sum(weights.*Ff);
    e       = Z(:,2)'*(weights.*Ff);
    sigmaE  = (c.*d-b.*e)./(a.*c-b.^2);
    sigmaA  = (-b.*d+a.*e)./(a.*c-b.^2);
    sigmaE(sigmaE<0)=0;sigmaA(sigmaA<0)=0;
    h2WLS   = sigmaA./(sigmaA+sigmaE);
    sigmaP  = sigmaA+sigmaE;
    h2WLS(Score<0)=0;
else
    h2WLS  = 0;
    sigmaA = 0;
    sigmaE = 0; 
end

varargout{1} = Ts;                  % Vectorized T_s
varargout{2} = sigmaA';              % Vectorized h2WLS  
varargout{3} = sigmaE';              % Vectorized sigmaP  
end
    