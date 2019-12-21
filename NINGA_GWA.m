function NINGA=NINGA_GWA(Y,M,i,Phi2,Opt)

if  ~isempty(Phi2)
    %Find the transformation matrix, apply it to gwa and data
    [vec,value]        = DimReduc(M,Phi2);
    %vec=csvread('vec.csv');
    NINGA.sy           = vec'*Y;size(NINGA.sy);
    %Auxiliary model design matrix
    Z                 = [ones(size(NINGA.sy,1),1) diag(value)]; 
    %calculate phenotypic variance
    sigmaP            =  mean(NINGA.sy.*NINGA.sy);
    %calculate random effects using 1-step estimator
    [~,sigmaA,sigmaE] = heterV(NINGA.sy,Z,1);
    sigmaE(sigmaA==0) = sigmaP(sigmaA==0);
    VC                = [sigmaE';sigmaA'];
    %find phenotypes that 1step performs bad and run fully convereged
    %estimator
    vv=(sigmaA+sigmaE)';Ind=find(vv>(1.5*sigmaP)| vv<(.5*sigmaP));
    if ~isempty(Ind)
        [~,vc,l]  = uni(NINGA.sy(:,Ind),value);
        VC(1,Ind) = vc(2,:);
        VC(2,Ind) = vc(1,:);
    end
    NINGA.h2          = VC(2,:)./sum(VC);
    NINGA.SigG        = VC(2,:);
    NINGA.SigE        = VC(1,:);
    weights           = 1./(Z*VC);
else
    [vec,value]        = DimReduc(M,[]);
    NINGA.sy           = vec'*Y;
    sigmaE             = mean(NINGA.sy.*NINGA.sy);
    weights            = repmat(1./sigmaE,size(NINGA.sy,1),1);
end
%load(fullfile('TMP',sprintf('Prep_chr_%s_%s_%s.mat',num2str(i),num2str(j),num2str(g))));
load(fullfile(Opt.GWAo,sprintf('Prep_chr_%s.mat',num2str(i))));

% %standardise the SNPs (make it compatible with FaST-LMM)
% nGc = size(gwa,2);nS=size(gwa,1);
% %gwa = ceil(gwa);
% for i=1:nGc
%     %calculate sample MAF
%     MAFC(i,1)    = (numel(find(gwa(:,i)==2))+.5*numel(find(gwa(:,i)==1)))/nS;
% end
% Pc   = repmat(MAFC',nS,1);
% gwa  = (gwa-2*Pc)./sqrt(2*Pc.*(1-Pc));
%Gmean  = repmat(mean(gwa),size(gwa,1),1);
%Gstd   = repmat(sqrt(var(gwa)),size(gwa,1),1);
%gwa    = (gwa-Gmean)./Gstd;
%clear Gmean Gstd
size(gwa)
sx              = single(vec'*gwa);
clear gwa
%TS vectors
tic
tmp             = single(weights.*NINGA.sy);
tmpX            = single((sx.^2)'); 
%NINGA.sy=[];
Score           = single(sx'*tmp);
Info            = single(tmpX*weights);
NINGA.TS        = single((Score.^2)./Info);
NINGA.BetaVar   = 1./Info;
NINGA.Beta      = (1./Info).*Score;
toc

tic
if ~isempty(Opt.clus)
    MaxC=zeros(size(TS,1),Opt.nP+1);
    for c=1:size(TS,1)
        [MaxC(c,1),NINGA.ClsSize{c,1},NINGA.ClsInd{c,1}]=NINGA_cluster(TS(c,:),Opt);
    end
end
toc
NINGA.MaxA(1,1) = max(NINGA.TS(:));

if Opt.nP>1
    if  ~isempty(Phi2)
         P           = csvread(fullfile(Opt.GWAo,'TMP_PermM'));
    else
         P          = csvread(fullfile(Opt.GWAo,'TMP_PermO'));
    end
    tic
    for p=1:Opt.nP
        Pvec                = P(:,p);
        %Pvec=randperm(size(NINGA.sy,1))';
        %size(Pvec)
        fprintf('(Chromosome %g: Permutation %g) \n',i,p)
        Info               = single(tmpX*weights(Pvec,:));
        Score              = single(sx'*tmp(Pvec,:));
        TS                 = single((Score.^2)./Info);
        NINGA.MaxA(p+1,1)   = max(TS(:));
        if ~isempty(Opt.clus)
            for c=1:size(TS,1)
                MaxC(c,p+1)=NINGA_cluster(TS(c,:),Opt);
            end
        end
    end
    if ~isempty(Opt.clus)
    NINGA.MaxC=max(MaxC)';
    end
    toc
end
end
