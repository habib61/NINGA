function tfce_img=NINGA_tfce(smoothed_img,E,H,Con,hmin,MaxE,type,Opt)
% Usage:
%tfce statisitc parameters should be specified based on image dimenstion.
%For instance for 2d images E=.5,H=1 and for 3d images H=2,E=.5
%
%Con defines connectivity criterion for clusetrs 
%Vanilla: tfce_img = tfce(smoothed_img,E,H,Con,hmin=0,hs=0,M1)
%M1     : tfce_img = tfce(smoothed_img,E,H,Con,hmin=0,hs=1.6,M1) 
%M3     : tfce_img = tfce(smoothed_img,E,H,Con,hmin=1.6,hs=0,M1) 
%M2     : tfce_img = tfce(smoothed_img,E,H,Con,hmin=0,hs=0,M2)
%M1_2   : tfce_img = tfce(smoothed_img,E,H,Con,hmin=0,hs=1.6,M2)
%M3_2   : tfce_img = tfce(smoothed_img,E,H,Con,hmin=1.6,hs=0,M2)
%M1a    :
%M3a    :
%
%_____________________________________
% Habib Ganjgahi
% Statistic Department. The univeristy of Warwick 
% January/2015


hmax     = max(smoothed_img(:)); 
deltah   = hmax/100;
tfce_img = zeros(size(smoothed_img));
       
switch type
    case 'M1'
        if  isempty(Opt.srf)
            %implement TFCE for volume
            for h = deltah:deltah:hmax
                %-Label each cluster with its own label using a connectivity criterion
                CC    = bwconncomp(smoothed_img>=h,Con);
                Size  = cellfun(@numel,CC.PixelIdxList);
                if ~isempty(MaxE)
                    M1      = min(Size,MaxE);
                else
                    M1      = Size;
                end
                integ = (M1.^E) * (h^H);
                for c = 1:CC.NumObjects
                    tfce_img(CC.PixelIdxList{c}) = ...
                        tfce_img(CC.PixelIdxList{c}) + integ(c);
                end
            end
        else
        %implement vanilla TFCE for surface
        % Vertexwise or facewise surface data
        tfce_img = zeros(size(smoothed_img));
        for h = deltah:deltah:max(smoothed_img(:))
            dpxl  = NINGA_dpxlabel(smoothed_img>=h,Opt.Yadjacency);
            U     = unique(dpxl(dpxl>0))';
            for u = 1:numel(U)
                idx = dpxl == U(u);
                tfce_img(idx) = tfce_img(idx) + ...
                    sum(Opt.Yarea(idx)).^Opt.tfce.E * h^Opt.tfce.H;
            end
        end
        end
        
        %TFCEM2
    case 'M2'
        for h = hmin:deltah:hmax
            %-Label each cluster with its own label using a connectivity criterion
            vol      = double(smoothed_img>=h);
            [c,num]  = spm_bwlabel(vol,Con);
           % Get size cluster size 
            if num>0
                clear ClusSize
               [ClusSize,ind]         = histc(c(:),(0:num) + 0.5);
               ClusSize               = ClusSize(1:end-1);
            %calculate e^H for all voxels at h level
                c(c>0)                = ClusSize(ind(ind>0));
                if ~isempty(MaxE)
                    M1                = min(c,MaxE);
                else
                    M1                = c;
                end
                tfce_img              = max(tfce_img,(h^H)*(M1.^E));
            end
        end
    case 'MS1'
        for h = hmin:deltah:hmax
            %-Label each cluster with its own label using a connectivity criterion           
            CC      = bwconncomp(smoothed_img>=h,Con);
           % Get size cluster size 
            if CC.NumObjects>0
               clear ClusSize
               ClusSize         = cellfun(@numel, CC.PixelIdxList);
               c                = labelmatrix(CC);
               ind              = c(cell2mat(CC.PixelIdxList'));
               ind2             = cell2mat(CC.PixelIdxList');
               c                = zeros(size(smoothed_img));
               c(ind2)          = double(ClusSize(ind));   
                if ~isempty(MaxE)
                    M1                = min(c,MaxE);
                else
                    M1                = c;
                end
                tfce_img              = tfce_img + ((h-hmin)^H)*(M1.^E);
            end
        end
    case 'MS2'
        for h = hmin:deltah:hmax
            %-Label each cluster with its own label using a connectivity criterion
            vol      = double(smoothed_img>=h);
            [c,num]  = spm_bwlabel(vol,Con);
           % Get size cluster size 
            if num>0
                clear ClusSize
               [ClusSize,ind]         = histc(c(:),(0:num) + 0.5);
               ClusSize               = ClusSize(1:end-1);
            %calculate e^H for all voxels at h level
                c(c>0)                = ClusSize(ind(ind>0));
                if ~isempty(MaxE)
                    M1                = min(c,MaxE);
                else
                    M1                = c;
                end
                tfce_img              = max(tfce_img,((h-hmin)^H)*(M1.^E));
            end
        end 
end
tfce_img = tfce_img *deltah;     
end