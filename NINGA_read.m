function [Y,readwith,extra,Mask] = NINGA_read(in,Mas)
%%This function reads input text file, recognizes type of data that can be
%%4D file or set of 3D images or a csv file. Also it checks for Mask 
%%availability, in case it doesn't exist it creats a mask. 
% 
%%
tmp              = textread(in,'%s');
[~,fname,fext]   = fileparts(tmp{1});
fext             = strdotsplit(strcat(fname,fext));
Mask.filename    = Mas;
%Read the mask if it's exist
%Check the input file
if numel(tmp)>1  %It means that the inputs are 3D image
    if strcmpi(fext{end-1},'nii') | strcmpi(fext{end},'gz')| strcmpi(fext{end},'nii')
        X        = NINGA_miscread(tmp{1},false);
        readwith = X.readwith;
        extra    = X.extra;
%         %Just4 Arash
%         tmpMask   = NINGA_miscread(Mask.filename,false);
%         Mask.data = tmpMask.data;
%         Y         = zeros(numel(tmp),numel(find(Mask.data>0)));
%         for r=1:numel(tmp)
%             X                = NINGA_miscread(tmp{r},false);
%             Y(r,:)           = X.data(Mask.data>0);
%         end
%        tmpData  = zeros([size(X.data) numel(tmp)]);
        %Read each 3D image and keep them as a 4D image to create Mask
        for r=1:numel(tmp)
            fprintf('(Reading the image from subject %g) \n',r)
            X                = NINGA_miscread(tmp{r},false);
            tmpData(:,:,:,r) = X.data;
        end    
        if ~isempty(Mas)
            tmpMask   = NINGA_miscread(Mask.filename,false);
            Mask.data = NINGA_MakeMask(tmpData).*tmpMask.data;
        else
            Mask.data = NINGA_MakeMask(tmpData);
        end
        %vectorize the image data
        Y        = zeros(numel(tmp),numel(find(Mask.data>0)));
        for r=1:numel(tmp)
            tt     = tmpData(:,:,:,r);
            Y(r,:) = tt(Mask.data>0);
        end
    %This means that the inputs are 3D images in HCP formats   
    else
        X        = NINGA_miscread(tmp{1},false,pwd);
        readwith = X.readwith;
        extra    = X.extra;
        if strcmp(readwith,'gifti')
            tmpData  = zeros([numel(tmp) size(X.data,2)]);
        else
            tmpData  = zeros([numel(tmp) size(X.data,1)]);
        end
        %Read each 3D image and vectorized it
        for r=1:numel(tmp)
            fprintf('(Reading the image from subject %g) \n',r)
            X            = NINGA_miscread(tmp{r},false,pwd);
            tmpData(r,:) = X.data;
        end

        Mask.data= NINGA_MakeMask(tmpData);
        Y        = tmpData(:,Mask.data>0);      
    end
elseif ~strcmp(fext,'csv') %it means that the input is either a 4D image or a csv file
    tmpData    = NINGA_miscread(tmp{1},false);
    readwith   = tmpData.readwith;
    extra      = tmpData.extra;
    %this means the input is HCP format, hence transpose
    if ndims(tmpData.data)==2
        tmpData.data = tmpData.data';    
    end
    if ~isempty(Mas)
        tmpMask     = NINGA_miscread(Mask.filename,false);
        Mask.data   = (tmpMask.data).*NINGA_MakeMask(tmpData.data);
    else
        Mask.data  = NINGA_MakeMask(tmpData.data);
    end
    %Vectorized the 4D data into a 2D matrix, rows are number of
    %subjects
    if ndims(tmpData.data)==2
        Y  = tmpData.data(:,Mask.data>0);
    else
        Y  = zeros(size(tmpData.data,4),length(find(Mask.data(:)>0)));
        for i=1:size(tmpData.data,4)
            kk              = tmpData.data(:,:,:,i);
            Y(i,:)          = kk(Mask.data>0);
        end
    end
else
    %the input is a csvfile
    Mask.data  = [];
    tmpData    = NINGA_miscread(tmp{1},false);
    Y          = tmpData.data;
    readwith   = tmpData.readwith;
    extra      = tmpData.extra;
end
    

    


% ==============================================================
function spl = strdotsplit(str)
% Split a string at the dots (.).
idx  = find(str == '.');
idxb = [1 idx+1];
idxe = [idx-1 numel(str)];
spl  = cell(numel(idxb),1);
for s = 1:numel(idxb),
    spl{s} = str(idxb(s):idxe(s));
end
