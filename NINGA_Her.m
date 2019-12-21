function NINGA = NINGA_Her(Y,M,S,Opt)

% 

%__________________________________________________________________________
%

% _________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% December/2014

%Output Initialisations
NINGA.Y.data     = cell(1,1);
NINGA.X          = cell(1,1);
NINGA.S          = cell(1,1);

NINGA.h2         = cell(1,1);
NINGA.h2Var      = cell(1,1);
NINGA.TS         = cell(1,1);
NINGA.SigmaY     = cell(1,1);


NINGA.Clus.Img   = cell(1,1);
NINGA.TFCE.Img   = cell(1,1);
NINGA.TFCE.FWE   = cell(1,1);
NINGA.Clus.thr   = Opt.clus;

NINGA.P.Uncor    = cell(1,1);
NINGA.P.Para     = cell(1,1);
NINGA.P.FWE      = cell(1,1);

switch Opt.Method
    case 'MLHer'
        %Calculate the EVD transformation matrix
        [vec,value]       = eig(S{1});
        sy                = vec'*Y;
        NINGA.sy          = sy;
        sx                = vec'*M;
        NINGA.sx          = sx;
        %Auxiliary model design matrix
        Z                 = [ones(size(sy,1),1) diag(value)]; 
        %Calculate the residuals
        nS                = size(sy,1);
        hat               = eye(nS)-sx*pinv(sx);  
        res               = hat*sy;
        beta              = pinv(sx)*sy;
        cov               = sx*beta;
    case 'REMLHer'
        [vec,value]       = DimReduc(M,S);
        res               = vec'*Y;
        sy                = res;
        NINGA.sy          = sy;
        sx                = [];
        Z                 = [ones(size(res,1),1) diag(value)]; 
end

%Calculate Unpermuted Score Statistic, and Sigma_A and SigmaE
[NINGA.TS{1},SigmaA,SigmaE]  = heterV(res,Z,1);

NINGA.h2{1}                  = SigmaA./(SigmaA+SigmaE);
NINGA.SigmaY{1}              = SigmaA+SigmaE;
NINGA.h2Var{1}               = NINGA_Var(Z,SigmaE,SigmaA)';

%It should be revise based on input type. Also based on Proper two
%sided permutation p-value
if Opt.Permout
    %Save permutation distribution
    if ~isempty(Opt.Mask)
        S          = palm_maskstruct(NINGA.TS{1},Opt.Mask.data,Opt.Mask.readwith,Opt.Mask.extra);
        S.filename = horzcat(Opt.o,'_',sprintf('p%.5d',0));
        palm_miscwrite(S)
    else
        fid   = fopen(sprintf('results_perm_%s.csv',Opt.o),'w');
        fprintf(fid,'c,');
        fprintf(fid,'%g,',NINGA.TS(:));
        fprintf(fid,'\n');
    end
end

if Opt.nP>0
    %Initialize output sturcture valraible
    %For Permutation
    IdS                 = zeros(size(NINGA.TS{1}));
    IdFWE               = zeros(size(NINGA.TS{1}));
    if Opt.v.I
        NINGA.v.MAX(1)  = max(NINGA.TS{1});
    end
    %Perform cluster-wise inference
    if ~isempty(NINGA.Clus.thr)  
        [NINGA.Clus.Max{1}(1),NINGA.Clus.Size{1}] = NINGA_cluster(NINGA.TS{1},Opt);
    end
    %Perform TFCE inference
    if Opt.tfce.I
        vol                   = zeros(size(Opt.Mask.data));
        %probability transform on chi images to perfom TFCE
        vol(Opt.Mask.data>0)  = norminv(chi2mixcdf(NINGA.TS{1},1));
        tmp                   = NINGA_tfce(vol,Opt.tfce.E,Opt.tfce.H,Opt.tfce.con,0,[],'M1',Opt);
        NINGA.TFCE.Img{1}     = tmp(Opt.Mask.data>0);
        NINGA.TFCE.Max(1)     = max(NINGA.TFCE.Img{1}(:));
        IdTFCEfwe             = zeros(size(NINGA.TFCE.Img{1}));
        IdTFCEun              = zeros(size(NINGA.TFCE.Img{1}));
    end
    for p=1:Opt.nP
        % Create permutation matrix for whole image
        %change the variable Perm to a structure and define different patterns
        fprintf('(Permutation %g) \n',p)
        switch Opt.Perm.method
            case 'Free'
                %Free Permutation
                [~,Pvec]    = sort(rand(size(res,1),1)); % CMC  
            case  'EX'
                %Restricted Permutation
                Opt.Perm.Xblk = FindEx(diag(value),[]);
                Pvec          = ExPerm(Opt.Perm.Xblk);
        end
        switch Opt.Method
            case 'MLHer'
                syP         = res(Pvec,:)+cov;
                resP        = hat*syP;
            case 'REMLHer'
                resP        = res(Pvec,:);
        end
        TSp         = heterV(resP,Z);
         %save the permutation distribution
        if Opt.Permout
            if ~isempty(Opt.Mask)
                S          = palm_maskstruct(TSp,Opt.Mask.data,Opt.Mask.readwith,Opt.Mask.extra);
                S.filename = horzcat(Opt.o,'_',sprintf('p%.5d',p));
                palm_miscwrite(S)
            else               
                fprintf(fid,'p%0.5d,',p);
                fprintf(fid,'%g,',TSp(:));
                fprintf(fid,'\n');
            end
        end
        %Indicator varaible to calculate uncorrected P-value
        IdS               = IdS+double(TSp>=NINGA.TS{1});
        IdFWE             = IdFWE+double((max(TSp)*ones(size(NINGA.TS{1})))>=NINGA.TS{1});
        if Opt.v.I
            NINGA.v.MAX(p+1,1)  = max(TSp);
        end   
        %Perform cluster-wise inference
        if ~isempty(NINGA.Clus.thr)
            NINGA.Clus.Max{1}(p+1,1) = NINGA_cluster(TSp,Opt);
        end
        %Perform TFCE inference
        if Opt.tfce.I
            vol                     = zeros(size(Opt.Mask.data));
            vol(Opt.Mask.data>0)    = norminv(chi2mixcdf(TSp,1));
            tmpTFCEp                = NINGA_tfce(vol,Opt.tfce.E,Opt.tfce.H,Opt.tfce.con,0,[],'M1',Opt);
            tmpTFCE                 = tmpTFCEp(Opt.Mask.data>0);
            NINGA.TFCE.Max(p+1,1)   = max(tmpTFCE(:));
            IdTFCEfwe               = IdTFCEfwe+double(max(tmpTFCE(:))*ones(size(NINGA.TFCE.Img{1}))>=NINGA.TFCE.Img{1});
            IdTFCEun                = IdTFCEun+double(tmpTFCE>=NINGA.TFCE.Img{1});
        end
    end
end
    if Opt.Permout && isempty(Opt.Mask)
        fclose(fid);
    end

    %Calculate Uncorretcted P-value
    if Opt.nP>0
        %calculate corrected and uncorrected voxel wise 1-P
        NINGA.P.Uncor{1}      = 1-((IdS+1)./(Opt.nP+1));
        NINGA.P.FWE{1}        = 1-(IdFWE+1)./(Opt.nP+1);
        %if requested, calculate corrected and uncorrected cluster size 1-P
        if ~isempty(Opt.clus)
             if ~isempty(Opt.clus)
                if isempty(NINGA.Clus.Size{1})
                    NINGA.Clus.FWEp     = 0;
                    NINGA.Clus.ImgP{1}  = zeros(size(Opt.Mask.data));
                else
                    for j=1:numel(NINGA.Clus.Size{1})
                        NINGA.Clus.FWEp{1}(1,j) = mean(NINGA.Clus.Size{1}(j)>=NINGA.Clus.Max{1});
                    end
                end      
             end 
        end
        %if requested, calculate corrected and uncorrected TFCE 1-P   
        if Opt.tfce.I
            NINGA.TFCE.FWE{1}   = zeros(size(Opt.Mask.data));
            NINGA.TFCE.uncor{1} = zeros(size(Opt.Mask.data));
            NINGA.TFCE.FWE{1}   = 1-(IdTFCEfwe+1)./(Opt.nP+1);
            NINGA.TFCE.uncor{1} = 1-((IdTFCEun+1)./(Opt.nP+1));
        end
    else
         NINGA.P.Uncor{1}     = [];
         NINGA.P.FWE{1}       = [];
    end
    NINGA.P.Para{1}     = chi2mixcdf(NINGA.TS{1},1);
    NINGA.Y.data{1}     = sy;
    NINGA.X             = sx; 
    NINGA.S             = diag(value);
end    