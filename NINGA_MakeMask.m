function mask = NINGA_MakeMask(tmpData)
% This function reads 4d or 2d images and make binary mask. 2d images
% correspond to HCP image formats.
% 
%
%
% _____________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% July/2015

%this is for removing constant voxels
removecte = true;

if ndims(tmpData)==4
    %make mask for 4d images 
    dim  = size(tmpData);
    mask = zeros(dim(1:3));
    for i=1:size(tmpData,2)
        for j=1:size(tmpData,3)
            kk          = squeeze(tmpData(:,i,j,:));
            IndNan      = any(isnan(kk),2);
            IndInf      = any(isinf(kk),2);
            IndCt       = sum(diff(kk,1,2).^2,2)==0;
            mask(:,i,j) = ~(IndNan | IndInf|IndCt);
        end
    end
else
    %make mask for 2d images
    ynan = any(isnan(tmpData),1);
    yinf = any(isinf(tmpData),1);
    if removecte,
        ycte = sum(diff(tmpData,1,1).^2) == 0;
    else
        ycte = false(size(yinf));
    end
    mask = ~ (ynan | yinf | ycte);
end
end