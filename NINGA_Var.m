function Var = NINGA_Var(Z,SigmaE,SigmaA)
% Calculates h2WLS variance

%__________________________________________________________________________
%

% _________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% December/2014



s       = [SigmaE' ;SigmaA'];
weights = 1./((Z*s).^2);
a       = sum(weights);
b       = Z(:,2)'*weights;
c       = (Z(:,2)'.^2)*weights;
det     = a.*c-b.^2;
G       = (SigmaA./((SigmaA+SigmaE).^2))';
E       = (SigmaE./((SigmaA+SigmaE).^2))';

Var     = 2*((G.^2).*c+2*(G.*E).*b+(E.^2).*a)./det;
end


