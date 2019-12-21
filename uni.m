function varargout = uni(varargin)
% This is a function to find the parameters
% Y: Data
% S: Correlation matrices for each component
% X: Design matrix with covariates (EVs)
%Habib Ganjgahi and Anderson Winkler


dflt.res = 5;
dflt.opt = optimset(        ...
    'Algorithm',    'SQP',  ...
    'Diagnostics',  'off',  ...
    'Display',      'off');


%call the actual variables not structures
fields = {'Y','val','X'};
for a = 1:nargin,
    eval(sprintf('%s = varargin{a};',fields{a}))
end

value   = diag(val);
[nS,nT] = size(Y);   % number of subjects and traits

%for the future find a way to this issue, number of VC should be found
%automaticly
%nK=1;
nK=2;

% % If an input design matrix has been input, use it. Otherwise
% % treat Y as already residualised
% if exist('X','var'),    
%     % Check if there is one X for all traits, of if X is trait-specific
%     if iscell(X),
%         Xone = false;
%         nC   = size(X{1},2);
%     else
%         Xone = true;
%         nC = size(X,2);
%     end
%     % Regress out using OLS to use as initial guesses for the MLE later
%     if Xone,
%         xtx  = diag(pinv(X'*X));
%         beta = X\Y;       % betas
%         r    = Y - X*beta;   % residuals
%     else
%         beta = cell(size(nT,1));
%         r    = zeros(size(Y));
%         xtx  = cell(size(nT,1));
%         for t = 1:nT,                          % for each trait
%             xtx{t}  = diag(pinv(X{t}'*X{t}));
%             beta{t} = X{t}\Y(:,t);             % betas
%             r(:,t)  = Y(:,t) - X{t}*beta{t};   % residuals
%         end
%     end
%     [p_wahh,sigma_wls,~] = heterR(Y,X,val,1,'wls');
% else
%     % If there are no design matrix, r = Y
%     r = Y;
%     nC = 0;
%     [p_wahh,sigma_wls,wahh] = heterR(Y,val,'wls');
% end
% % Find the initial value via WLS heritability estimator
% initK                   = [p_wahh;sigma_wls];
nC=0;
[~,sigmaA,sigmaE] = heterV(Y,[ones(nS,1) value],1);
initK             = [sigmaA';sigmaE'];
initK =repmat([.5 ;.5],1,nT);
fprintf(1,'Optimising parameters: [ ');
% Linear constraits and optimisation options
%Aeq = [zeros(1,nC) ones(1,nK)];  % equality constraints only
%beq = 1;
% Init some vars for the for-loop below
parms = zeros(nC+nK,nT); % betas and diagonal elements of R (univariate)
LL1   = zeros(1,nT);     % loglikelihoods for the correct model


% 
for t = 1:nT,
    fprintf(1,'%d',t);
    if nC > 0,
        % one design matrix per trait
        lb   = [beta{t}-(xtx{t}*40); zeros(nK,1)];
        ub   = [beta{t}+(xtx{t}*40); ones(nK,1)*Inf];
        init = [beta{t}; initK(:,t)];
        noncl= @(x) deal(myconstr(x),[]);
        objH = @(x)logi(Y(:,t),X{t},x,value);
     else % no design matrix supplied
         lb   = zeros(nK,1);
         ub   = ones(nK,1);
         %ub   = ones(nK,1)*Inf;
         init = initK(:,t);
         noncl= @(x)deal(myconstr2(x),[]);
         objH = @(x)logi(Y(:,t),x,value); 
    end
    fprintf(1,' ');  
    % Optimisation
    if isoctave,
        eqH = @(x)equality(Aeq,beq,x); % equality constraint handle
        [parms(:,t),LL1(t)] = sqp( ...
            init,         ...  % initialization values
            objH,         ...  % objective function handle
            eqH,          ...  % equality handle
            [],           ...  % inequality handle
            lb,ub,        ...  % lower & upper bounds
            dflt.opt.MaxIter); % options
    else
        [parms(:,t),LL1(t)] = fmincon( ...
            objH,         ...  % objective function handle
            init,         ...  % initialization values
            [],[],     ...  % inequality constraints
            [],[],      ...  % equality constraints
            lb,ub,        ...  % lower & upper bounds
            noncl,dflt.opt);      % options
    end
end
fprintf(1,']\n');
% Output results

varargout{1} = parms(1:nC,:);                % betas (MLE)          
varargout{2} = parms(nC+1:nC+nK,:);          % varcomps
varargout{3} = LL1;                          % loglikes
