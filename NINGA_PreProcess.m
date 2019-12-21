function  [Y,ped,M,C,Opt]=NINGA_PreProcess(Opt)

%% 
%% Read the data
[Y,readwith,extra,Opt.Mask] = NINGA_read(Opt.in,Opt.Mask.filename);
Opt.Mask.readwith           = readwith;
Opt.Mask.extra              = extra;

%% Read the design matrix and contrast file 
[~,fname,fext]   = fileparts(Opt.d);
if strcmp(fext,'.mat')
    M = NINGA_vestread(Opt.d);
else
    M = csvread(Opt.d);
end

%read the contrast file if exist 
if isempty(Opt.c)
    Opt.nCon = 1;
    C        = [];
else
    [~,fname,fext] = fileparts(Opt.c);
    if strcmp(fext,'.con')
        C          = NINGA_vestread(Opt.c);
    else
        C          = csvread(Opt.c);
    end
    Opt.nCon       = size(C,1);
end

%% Read the kinship and covert it to a square matrix
ped                         = NINGA_kin(Opt.ped,Opt.kin);


%% read surface file for cluster or tfce statistics
%at the moment it works only for vertexdata not face

if ~isempty(Opt.srf)
    Opt.Yisvtx     = true;
    Opt.srf        = NINGA_miscread(Opt.srf);
    Opt.Yadjacency = NINGA_adjacency(Opt.srf.data.fac,Opt.Yisvtx);
    % Load areas
    if isempty(Opt.srfare)
       % No area means area of the actual surface file
       Opt.srfarea.data = [];
    elseif exist(Opt.srfare,'file')
        % A file with the average areas from native geometry
        Opt.srfarea = NINGA_miscread(Opt.srfare,false);
    elseif ~ isnan(str2double(Opt.srfare))
        % A weight (such as 1)
        Opt.srfarea.data = str2double(Opt.srfare);
    else
       error('Unknown option given to "Opt.srf" or file doesn''t exist:\n'); 
    end
        % Surface area, to be used by the spatial statistics
    if isempty(Opt.srfarea.data),
        % No area file given, use the actual surface area
        Opt.Yarea  = NINGA_calcarea(Opt.srf.data,Opt.Yisvtx);
    elseif numel(Opt.srfarea.data) == 1,
        % A weight given (such as 1): use that for each vertex or face,
        % treating as if all had the same area.
        if Opt.Yisvtx,
            Opt.Yarea = Opt.srfarea.data .* ...
                ones(size(Opt.srf.data.vtx,1),1);
        elseif Opt.Yisfac,
            Opt.Yarea = Opt.srfarea.data .* ...
                ones(size(Opt.srf.data.fac,1),1);
        end

    else
        % Otherwise, just use the data from the file (already loaded).
        Opt.Yarea = Opt.srfarea.data;
    end
end

%% Read subject ID and convert it to char
SubjID                      = strcsvread(Opt.subjID);

%now remove extra subjects and reshpae phenotype based on kinship id
if ~strcmpi(Opt.Method,'GWAS')
    %find subjects that doesn't have phenotype
    KinId      = cellfun(@num2str,ped.id,'UniformOutput',false);
    SubjID     = cellfun(@num2str,SubjID,'UniformOutput',false);
    delDataKin = setxor(SubjID ,KinId);
    idxKin     = false(size(KinId));
    idxData2   = false(size(SubjID));
    for s = 1:numel(delDataKin)
        idxKin   = idxKin  | strcmpi(delDataKin{s},KinId);
        idxData2 = idxData2| strcmpi(delDataKin{s},SubjID);
    end
    ped.id(idxKin)    = []; 
    KinId(idxKin)     = [];
    ped.phi2          = ped.phi2(~idxKin,~idxKin);
    
    SubjID(idxData2)  = [];
    Y(idxData2,:)     = [];
    M(idxData2,:)     = [];
   
    %find rearranging indices with respect to kinID
    nS      = numel(KinId);
    IndKin  = zeros(nS,1);
    IndData = zeros(nS,1);
    for l=1:nS
        IndData(l)  = find(cell2mat(cellfun(@(x) strcmpi(x,KinId{l}), SubjID,'UniformOutput',false)));
    end
    fprintf('Matching the individuals between the kinship id and Subject id and Rearranging them based kinship id \n')
    Y          = Y(IndData,:);
    M          = M(IndData,:);
    SubjID     = SubjID(IndData,:);
   
    if strcmp(readwith,'load')
       Opt.Mask.data=[];
    end
%     fprintf('Saving the Rearranged data \n')
%     if strcmp(readwith,'wb_command') | strcmp(readwith,'gifti')
%         Y=Y';
%     end
%     S          = NINGA_maskstruct(Y,Opt.Mask.data,readwith,extra);
%     S.filename = 'Prep_Data';
%     NINGA_miscwrite(S)
%     
%     csvwrite('Prep_design',M);
%     cellcsvwrite({fullfile(pwd,strcat('Prep_Data.',Opt.Mask.extra.cifti_file_extension,'.nii'))},'Prep_Data.txt')
%     cellcsvwrite({fullfile(pwd,strcat('Prep_Data.',Opt.Mask.readwith,'.nii'))},'Prep_Data.txt')
%     csvwrite('subjid_reorder.csv',SubjID);
%     save('Prep_data','ped','M')   
    %save('Prep_data','ped','M','Y','SubjID','readwith','extra')   
else
    %Read GWA ID
    GWAid                      = strcsvread(Opt.GWAid);
    KinID                      = ped.id;
    
    %Find the subjects that have all kin,data and gwa info and
    %rearrange them based in gwa order
    fprintf('Pre-processing: Intersection ....')
    [IndKin,IndData,idxData,idxGWA,idxKin,idxData2] = ...
                  NINGA_Reshape(SubjID,KinID,GWAid);
    fprintf(' done \n')                         
    %remove extra subjects and rearrange them based on GWAid  
    fprintf('removing extra subjects and rearrange them based on GWAid') 
    SubjID(idxData)   = [];
    Y(idxData,:)      = [];
    M(idxData,:)      = [];
    SubjID(idxData2)  = [];       
    Y(idxData2,:)     = [];
    M(idxData2,:)     = [];

    %remove extra subjects from kinship and rearrange them
    KinID(idxKin)     = [];
    ped.id(idxKin)    = [];
    ped.phi2          = ped.phi2(~idxKin,~idxKin);
    ped.phi2          = ped.phi2(IndKin,IndKin);

    Y                 = Y(IndData,:);
    M                 = M(IndData,:);
    SubjID            = SubjID(IndData); 
    if ~exist(Opt.GWAo,'file')
         cmd = sprintf('mkdir %s',Opt.GWAo);
         system(cmd);
    end
    i                 = str2num(getenv('SGE_TASK_ID'));
    %i=1;
     %Generate the permutation matrix
    if i==1
        if Opt.nP>0
            if Opt.inormal
                [vec,value] = DimReduc(M(:,2:end),ped.phi2);
            else
                [vec,value] = DimReduc(M,ped.phi2);
            end
            nS              = size(vec',1);
            for p=1:Opt.nP 
                Pvec(:,p)   = randperm(nS)';
            end
            csvwrite(fullfile(Opt.GWAo,'TMP_PermM'),Pvec)
        end
    end 
    
     %remove extra subjects from gwa
     
     %gwa=load(fullfile(Opt.GWAin,sprintf('%s_%s.mat',Opt.GWAb,num2str(i))));
     load(fullfile(Opt.GWAin,sprintf('%s_%s.mat',Opt.GWAb,num2str(i))));
     gwa=dosage;
     clear dosage
     gwa(idxGWA,:)     = [];
     gwa(idxData2,:)   = [];
     %save the mathced gwa 
     save(fullfile(Opt.GWAo,sprintf('Prep_chr_%s.mat',num2str(i))),'gwa','-v7.3');

    Phi2 = ped.phi2; Mask=Opt.Mask;
    save(fullfile(Opt.GWAo,sprintf('Misc_%s.mat',num2str(i))),'Phi2','Y','M','Mask','-v7.3');
end
end
