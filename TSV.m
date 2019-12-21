function  varargout=TSV(varargin)
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

fields = {'sy','sx','weights'};
for a = 1:nargin,
    eval(sprintf('%s = varargin{a};',fields{a}))
end


%calculate the score function
score  =  sx'*(weights.*sy);

%Calculate information matrix
infor  = (sx.*sx)'*weights;

%Calculate Score test 
TS     = (score.*score)./infor;
       
varargout{1} = TS;                % Vectorized T_s
end