function  varargout=heterR(varargin)
%%check for the covariate matrix 
%echo on
if nargin == 3 || nargin == 4
    fields = {'sy','value','test','Pvec'};
    for a = 1:nargin,
    eval(sprintf('%s = varargin{a};',fields{a}))
    end
    %check if design matrix is provided
    nT = size(sy,2);
    for t=1:nT
        f(:,t)=sy(:,t).^2; 
    end
else nargin == 5 || nargin == 6;
     %if the covarite matrix is provided
     fields = {'sy','sx','value','k0','test','Pvec'};
        for a = 1:nargin,
           eval(sprintf('%s = varargin{a};',fields{a}))
        end
        nT = size(sy,2);
    for t=1:nT
        beta_s{t} = sx{t}\sy(:,t);               % Calculate residuals
        rs(:,t)   = (sy(:,t) - sx{t}*beta_s{t});   
        f(:,t)    = rs(:,t).^2;                  % Calculate Squared residuals
    end
end 
        
[nS,nT] = size(f);


z=diag(value);
%permute the eigenvalues
if exist('Pvec','var')
    z=z(Pvec);
end
z1      = [ones(nS,1) z];

for t=1:nT
    sigma     = mean(f(:,t));                % Esitmate phenotypic variance
    switch test;
        case 'bp'
            %sprintf('Calculating the score test for trait %s',num2str(t))  
             %Calculate the score function
             if 1/(2*sigma)*z'*(f(:,t)./sigma-ones(nS,1))<=0
                    s(:,t)  = [0;0];
                    h(:,t)  = 0;
                    LM(:,t) = 0;
                    sig(:,t)= sigma;
             else
                    % Auxiliary model (regress squared residuals on Eigenvalues minus 1)                   
                    s(:,t)= z1\f(:,t);
                    %use the constaint least square if the intercept or
                    %slope of the auxiliary model is negative
                    if s(1,t)<0 | s(2,t)<0
                        s(:,t)=lsqnonneg(z1,f(:,t));
                    end
                    h(:,t)       = s(2,t);
                    sig(:,t)     = s(1,t);
                    sigma_inv    = 1/sigma;
                    %calcute the score statistics
                    z0        = sigma_inv*z;
                    D         = eye(nS)-sigma_inv*ones(nS,1)*inv(ones(1,nS)*sigma_inv.^2*ones(nS,1))*ones(1,nS)*sigma_inv;
                    LM(:,t)   = .5*s(2,t)'*(z0'*D*z0)*s(2,t);    %score test  
             end
         case 'wls';
            %sprintf('Calculating the Wald test for trait %s in term of WLS estimator' ,num2str(t))
             %check the score function, if its value at zero is
             %negative then WLS estimator and the wald test is zero
             if 1/(2*sigma)*z'*(f(:,t)./sigma-ones(nS,1))<=0
                 alphaz(:,t) = [0;0];
                 sig(:,t)    = mean(f(:,t)); 
                 h(:,t)      = 0;
                 LM(:,t)     = 0;
             else
                 %calculate the weights in terms of OLS
                 s(:,t)= z1\f(:,t);
                 if s(1,t)<0 | s(2,t)<0
                    s(:,t)=lsqnonneg(z1,f(:,t));    
                 end
                 %Use WLS to estimate the auxiliary model parameters
                 sigma_inv   = diag(1./(z1*s(:,t))); 
                 alphaz(:,t) = inv(z1'*sigma_inv.^2*z1)*z1'*(sigma_inv.^2)*f(:,t);   
                 %alphaz(:,t) = pinv(z1'*sigma_inv.^2*z1)*z1'*(sigma_inv.^2)*f(:,t); 
                 %Check the estimated parameters and constrain the wrong
                 %values
                 if alphaz(1,t)<0
                    alphaz(1,t)=0;
                 end
                 if alphaz(2,t)<0
                    alphaz(2,t)=0;
                 end
                 h(:,t)      = alphaz(2,t);
                 sig(:,t)    = alphaz(1,t); 
              
                 %Calculate the Wald test in terms of WLS estimator
                 %consdier inv below to remove it
                 z0      = sigma_inv*z;
                 D       = eye(nS)-sigma_inv*ones(nS,1)*inv(ones(1,nS)*sigma_inv.^2*ones(nS,1))*ones(1,nS)*sigma_inv;
                 LM(:,t) = 0.5*alphaz(2,t)'*(z0'*D*z0)*alphaz(2,t);
             end
        case 'wols';
            if 1/(2*sigma)*z'*(f(:,t)./sigma-ones(nS,1))<=0
                h(:,t)   = 0;
                sig(:,t)  = 0;
                LM(:,t)   = 0;
            else
                s(:,t)= z1\f(:,t);
                 if s(1,t)<0 | s(2,t)<0
                    s(:,t)=lsqnonneg(z1,f(:,t));    
                 end
               sigma_invs = 1/mean(f(:,t));
               %calcute the score statistics
               z0s        = sigma_invs*z;
               Ds         = eye(nS)-sigma_invs*ones(nS,1)*inv(ones(1,nS)*sigma_invs.^2*ones(nS,1))*ones(1,nS)*sigma_invs;
               LM(:,t)    = .5*s(2,t)'*(z0s'*Ds*z0s)*s(2,t);    %score test
               h(:,t)     = LM(:,t);
               %calculate the wls
               %covariance matrix of the first iteration
               sigma_inv  = diag(1./(z1*s(:,t)));
               % WLS paramet estimation
               alphaz(:,t) = inv(z1'*sigma_inv.^2*z1)*z1'*(sigma_inv.^2)*f(:,t);
               if alphaz(1,t)<0
                    alphaz(1,t)=0;
                 end
                 if alphaz(2,t)<0
                    alphaz(2,t)=0;
                 end
               %Calculate the Wald test in terms of WLS estimator
               %consdier inv below to remove it
               z0         = sigma_inv*z;
               D          = eye(nS)-sigma_inv*ones(nS,1)*inv(ones(1,nS)*sigma_inv.^2*ones(nS,1))*ones(1,nS)*sigma_inv;
               wahh(:,t)  = .5*alphaz(2,t)'*(z0'*D*z0)*alphaz(2,t);
               sig(:,t)    = wahh(:,t);      
            end
    end
end
varargout{1} = h;                % estimated h2
varargout{2} = sig;              % estimated phenotypic variance  
varargout{3} = LM;               % test statisitcs
%LM=1-chi2cdf(Lm,1);
end