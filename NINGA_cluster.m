function [Max,ClsSize]=NINGA_cluster(TS,Opt)
%% Returns maximum statistic and maximum cluster size for a given cluster
% forming threshold
%
% Usage [max_s,MaxClVox,ClSize] = max_clsf(TS,e,Dim,thresh);
% - TS      - Test statistic image
% - Mask    - voxel indixies in image space
% - Dim     - statistic image dimension
% - thr     - Cluster forming threshold in terms of P-value like P=0.01
%
% Outputs
% max_s     - Maximum statistic
% MaxClVox  - Maximum cluster size for a given cluster forming threshold
% ClSize    - cluster size distribution for a given threshold
%__________________________________________________________________________
%
% max_clsf find maximum statistic value in the statistic image. Also it
% finds clusters, their sizes and returns maximum cluster size for the
% given cluster forming threshold.
% 
%__________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% December/2014
%%

if ndims(TS)>2
    vol = TS;
else
    vol                   = zeros(size(Opt.Mask.data));
    vol(Opt.Mask.data>0)  = TS;
end

if isempty(Opt.srf)
    %Find clusters and calculate their size, not bwconncomp is faster 
    CC      = bwconncomp(vol>=Opt.clus,Opt.con);
    if CC.NumObjects>0
        Size     = cellfun(@numel, CC.PixelIdxList);
        Max      = max(Size);
        if nargout > 1
           ClsSize  = zeros(size(vol));
           for c = 1:CC.NumObjects
                ClsSize(CC.PixelIdxList{c}) = Size(c);
           end
           %vectorized cluster size image
           ClsSize  = ClsSize(Opt.Mask.data>0)';
        end
    else
       ClsSize = [];
       Max     = 0;
    end
else
    Dt    = vol>=Opt.clus;
    if sum(Dt>0)
        % Connected components:
        dpxl  = NINGA_dpxlabel(Dt,Opt.Yadjacency);
        % Compute the cluster stats
        U     = unique(dpxl(dpxl>0))';
        Size  = zeros(size(U));
        for u = 1:numel(U)
            Size(u) = sum(Opt.Yarea(dpxl == U(u)));
        end
        Max     = max(Size);
        % Compute the statistic image (this is normally for the 1st perm only)
        if nargout > 1
            ClsSize = zeros(size(Dt));
            for u = 1:numel(U)
                ClsSize(dpxl == U(u)) = Size(u);
            end
            ClsSize = ClsSize(Opt.Mask.data>0)';
        end
    else
        ClsSize = [];
        Max     = 0;
    end
end
    
end

