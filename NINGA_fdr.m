function padj = NINGA_fdr(pval)
% Compute FDR-adjusted p-values

V           = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR]   = sort(oidx);
padj        = zeros(size(pval));
prev        = 1;
for i = V:-1:1,
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);