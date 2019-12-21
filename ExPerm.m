function Pvec = ExPerm(Xblk)
%Number of Subjects
nScan = sum(Xblk);
%Naive basic idea to find exchangeable blocks
% %size of exchangeable blocks
% %Xblk  = [4 3 5]
% kk  = floor(10*diag(value));
% Ind = unique(kk);

CsXblk = cumsum(Xblk);
%Number of exchangeable blocks
nXblk    = length(Xblk);

%-Generate a new random permutation - see randperm
p = cell(1,nXblk);
for i=1:nXblk
    [~,p{i}] = sort(rand(1,Xblk(i))); 
    if i>1
        p{i}     = p{i} + CsXblk(i-1);
    end
end
Pvec = (cell2mat(p))';
end
