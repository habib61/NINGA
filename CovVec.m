function NINGA=CovVec(Y,M,C,Ss,Opt)
% 

%__________________________________________________________________________
%

% _________________________________________________________________________
% Habib Ganjgahi
% Statistic Department, The univeristy of Warwick 
% December/2014

%Output Initialisations
NINGA.Y.data     = cell(Opt.nCon,1);
NINGA.X          = cell(Opt.nCon,1);
NINGA.S          = cell(Opt.nCon,1);

NINGA.Beta       = cell(Opt.nCon,1);
NINGA.BetaVar    = cell(Opt.nCon,1);
NINGA.TS         = cell(Opt.nCon,1);
NINGA.SigmaY     = cell(Opt.nCon,1);

NINGA.UncorP     = cell(Opt.nCon,1);
NINGA.MAXv       = cell(Opt.nCon,1);
NINGA.Clus.Img   = cell(Opt.nCon,1);
NINGA.TFCE.Img   = cell(Opt.nCon,1);



NINGA.Clus.thr   = Opt.clus;
for c=1:Opt.nCon
    %Ortogonalize the design with respecto to the contrast
    [X0,Zz,eCm,eCx]       = NINGA_partition(M,C(c,:)');
    %unpermutaed data
    %corraleted data 
    if  ~isempty(Ss)
        %In case of related subjects, find the projection matrix
        [vec,value]       = DimReduc(Zz,Ss);
        sy                = vec'*Y;
        sx                = vec'*X0;
        %Auxiliary model design matrix
        Z                 = [ones(size(sy,1),1) diag(value)]; 
        %TS vectors
        %calculate \Sigma^-1 by 1-step estimator
        [~,sigmaA,sigmaE] = heterV(sy,Z,1);
        VC                = [sigmaE';sigmaA'];
        NINGA.h2          = sigmaA./(sigmaA+sigmaE);
        NINGA.sigmaA      = sigmaA;NINGA.sigmaE=sigmaE;
        %calculate weights
        weights           = 1./(Z*VC);
        %Beta vector
        NINGA.Beta{c}      = BV(sy,sx,Z);
        NINGA.BetaVar{c}   = VarB(sy,sx,Z);
        %Estimate ReML phenotype Covaraince
        hatF              = eye(size(sx,1))-sx*inv(sx'*sx)*sx';
        hatF              = (hatF+hatF')/2;
        syF               = hatF*sy;
        [~,sigmaA,sigmaE] = heterV(syF,Z,1);
        VCF               = [sigmaE';sigmaA'];
        %phenotypic varaince
        %NINGA.SigmaY{c}    = Z*VCF;
    else
    
        %In case of unrelated subjects
        [vec,value]       = DimReduc(Zz,[]);
        sy                = vec'*Y;
        sx                = vec'*X0;
        NINGA.Beta{c}      = sx\sy;
        NINGA.BetaVar{c}   = VarB(sy,sx,[]);
        sigmaE            = mean(sy.*sy);
        weights           = repmat(1./sigmaE,size(sy,1),1);
        %Estimate ReML phenotype Covaraince
        hatF              = eye(size(sx,1))-sx*inv(sx'*sx)*sx';
        hatF              = (hatF+hatF')/2;
        syF               = hatF*sy;
        %phenotypic varaince
        %NINGA.SigmaY{c}   = repmat(mean(syF.*syF),size(sy,1),1);
    end

    %Calculate vectorized score test
    Score                 = sx'*(weights.*sy);
    NINGA.TS{c}           = sign(Score).*sqrt(TSV(sy,sx,weights));
    
    %It should be revise based on input type. Alos based on Proper two
    %sided permutation p-value
    if Opt.Permout
        %Save permutation distribution
        if ~isempty(Opt.Mask)
            S          = palm_maskstruct(NINGA.TS{c},Opt.Mask.data,Opt.Mask.readwith,Opt.Mask.extra);
            S.filename = horzcat(Opt.o,'_',num2str(c),'_',sprintf('p%.5d',0));
            palm_miscwrite(S)
        else
            fid   = fopen(sprintf('results_perm_%s.csv',Opt.o),'w');
            fprintf(fid,'c,');
            fprintf(fid,'%g,',NINGA.TS{c}(:));
            fprintf(fid,'\n');
        end
    end
    if Opt.nP>0
        %Initialize output sturcture valraible
        %For Permutation
        IdS                 = zeros(size(NINGA.TS{c}));
        IdFWE               = zeros(size(NINGA.TS{c}));
        %rng(Opt.seed)
        if Opt.v.I
            NINGA.v.MAX{c}(1)  = max(NINGA.TS{c});
        end
        %Perform cluster-wise inference
        if ~isempty(NINGA.Clus.thr)
           [NINGA.Clus.Max{c}(1),NINGA.Clus.Size{c}] = NINGA_cluster(NINGA.TS{c},Opt);
        end
        %Perform TFCE inference
        tic
        if Opt.tfce.I
            vol                   = zeros(size(Opt.Mask.data));
            vol(Opt.Mask.data>0)  = NINGA.TS{c};
            tmp                   = NINGA_tfce(vol,Opt.tfce.E,Opt.tfce.H,Opt.tfce.con,0,[],'M1',Opt);
            NINGA.TFCE.Img{c}     = tmp(Opt.Mask.data>0);
            NINGA.TFCE.Max{c}(1)  = max(NINGA.TFCE.Img{c}(:));
            IdTFCEfwe             = zeros(size(NINGA.TFCE.Img{c}));
            IdTFCEun              = zeros(size(NINGA.TFCE.Img{c}));
        end
        toc
        for p=1:Opt.nP
            % Create permutation matrix for whole image
            %change the variable Perm to a structure and define different patterns
            fprintf('(Permutation %g) \n',p)
            switch Opt.Perm.method
                case 'Free'
                    %Free Permutation
                    [~,Pvec]    = sort(rand(size(sy,1),1)); % CMC  
                case  'EX'
                    %Restricted Permutation
                    Opt.Perm.Xblk = FindEx(diag(value),[]);
                    Pvec          = ExPerm(Opt.Perm.Xblk);
            end
            syP         = sy(Pvec,:);
            ScoreP      = sx'*(weights(Pvec,:).*syP);
            TSp         = sign(ScoreP).*sqrt(TSV(syP,sx,weights(Pvec,:)));
             %save the permutation distribution
            if Opt.Permout
                if ~isempty(Opt.Mask)
                    S          = palm_maskstruct(TSp,Opt.Mask.data,Opt.Mask.readwith,Opt.Mask.extra);
                    S.filename = horzcat(Opt.o,'_',num2str(c),'_',sprintf('p%.5d',p));
                    palm_miscwrite(S)
                else               
                    fprintf(fid,'p%0.5d,',p);
                    fprintf(fid,'%g,',TSp(:));
                    fprintf(fid,'\n');
                end
            end
            
            %Indicator varaible to calculate uncorrected P-value
            IdS               = IdS+double(TSp>=NINGA.TS{c});
            IdFWE             = IdFWE+double((max(TSp)*ones(size(NINGA.TS{c})))>=NINGA.TS{c});
            if Opt.v.I
                NINGA.v.MAX{c}(p+1,1)  = max(TSp);
            end
            
            %Perform cluster-wise inference
            if ~isempty(NINGA.Clus.thr)
                NINGA.Clus.Max{c}(p+1,1) = NINGA_cluster(TSp,Opt);
            end
            %Perform TFCE inference
            if Opt.tfce.I
                vol                      = zeros(size(Opt.Mask.data));
                vol(Opt.Mask.data>0)     = TSp;
                tmpTFCEp                 = NINGA_tfce(vol,Opt.tfce.E,Opt.tfce.H,Opt.tfce.con,0,[],'M1',Opt);
                tmpTFCE                  = tmpTFCEp(Opt.Mask.data>0);
                NINGA.TFCE.Max{c}(p+1,1) = max(tmpTFCE(:));
                IdTFCEfwe                = IdTFCEfwe+double(max(tmpTFCE(:))*ones(size(NINGA.TFCE.Img{c}))>=NINGA.TFCE.Img{c});
                IdTFCEun                 = IdTFCEun+double(tmpTFCE>=NINGA.TFCE.Img{c});
            end
        end
        %Calculate Uncorretcted P-value
        if Opt.nP>0
            %calculate corrected and uncorrected voxel-wise 1 - P-values
            NINGA.P.Uncor{c}  = 1-((IdS+1)./(Opt.nP+1));
            NINGA.P.FWE{c}    = 1-(IdFWE+1)./(Opt.nP+1);
            %if requested, calculate cluster size FWE 1 - P-values
            if ~isempty(Opt.clus)
                 if ~isempty(Opt.clus)
                    if isempty(NINGA.Clus.Size{c})
                        NINGA.Clus.FWEp     = 0;
                        NINGA.Clus.ImgP{c}  = zeros(size(Opt.Mask.data));
                    else
                        for j=1:numel(NINGA.Clus.Size{c})
                            NINGA.Clus.FWEp{c}(1,j) = mean(NINGA.Clus.Size{c}(j)>=NINGA.Clus.Max{c});
                        end
                    end      
                end     
            end
            % if requested calculated corrected and uncorrected TFCE 1 - P-values 
            if Opt.tfce.I
                NINGA.TFCE.FWE{c}   = zeros(size(Opt.Mask.data));
                NINGA.TFCE.uncor{c} = zeros(size(Opt.Mask.data));
                NINGA.TFCE.FWE{c}   = 1-(IdTFCEfwe+1)./(Opt.nP+1);
                NINGA.TFCE.uncor{c} = 1-((IdTFCEun+1)./(Opt.nP+1));
            end
        else
            NINGA.P.Uncor{c}  = [];
            NINGA.P.FWE{c}    = [];
        end
    end
    if Opt.Permout && isempty(Opt.Mask)
        fclose(fid);
    end
    NINGA.P.Para{c}  = normcdf(NINGA.TS{c});
    NINGA.X{c}       = sx;
    NINGA.S{c}       = diag(value);
end

end
