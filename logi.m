function varargout = logi(varargin)
if nargin == 3,
    % Call "L = loglike(Y,bb,value)"
    % Accept arguments
    fields = {'Y','bb','value'};
    for a = 1:nargin,
        v.(fields{a}) = varargin{a};
    end
    n=size(v.value,1);
    nC=0;
    % Residulas
    tau = v.Y;
    c   = v.bb(nC+1)*(v.value)+v.bb(nC+2)*ones(n,1);
    %gradL(1)=.5*(sum((v.value-1)./c)+(1/v.f1^2)*sum((tau.^2.*(-v.value+1))./c));
elseif nargin == 4,
    % Call "L = loglike(Y,X,bb,value)"
    % Accept arguments
    fields = {'Y','X','bb','value'};
    for a = 1:nargin,
        v.(fields{a}) = varargin{a};
    end
    % Make the residuals
    nC = size(v.X,2);
    n  = size(v.Y,1);
    %tau = (v.Y - v.X * v.bb(1:nC))./v.f1;
    tau = (v.Y - v.X * v.bb(1:nC));
    c   = v.bb(nC+1)*(v.value)+v.bb(nC+2)*ones(n,1);
%     gradL(1)=sum((tau.*v.X)./c)*.5*v.f1^2;
%     gradL(2)=.5*(sum((v.value-1)./c)+(1/v.f1^2)*sum((tau.^2.*(-v.value+1))./c));
end
 
 b=sum(log(c));
 a=sum(tau.^2./c);

 L=(a+b+numel(tau)*log(2*pi))/2;
 varargout{1}=L;
 %varargout{2}=gradL;
 
