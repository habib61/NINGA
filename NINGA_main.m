function [NINGA,Opt] = NINGA_main
% set up default parameters

%
%
% _____________________________________
% Habib Ganjgahi
% Statistic Department, uni of Oxford 
% Dec/2016
%%
Opt             = NINGA_defaults;

%% read the data
if ~Opt.PreProcess
    [Y,ped,M,C,Opt]             = NINGA_PreProcess(Opt);
else
    [Y,readwith,extra,Opt.Mask] = NINGA_read('Prep_Data.txt',Opt.Mask.filename);
    Opt.Mask.readwith           = readwith;
    Opt.Mask.extra              = extra;
    load('Prep_data')      
end

if Opt.inormal
     Y       = inormal(Y,'solar');
%    Y       = Y-repmat(mean(Y),size(Y,1),1);
%    M       = M-repmat(mean(M),size(M,1),1);
%    M (:,1) = [];
%    b       = M\Y;
%    cov     = M*b;
%    Y       = inormal(Y-cov,'solar')+cov;
end
%% Actual analysis starts here
switch Opt.Method
    case 'MLHer'
        %Call NINGA_Her to run ML heritability analysis
        NINGA          = NINGA_Her(Y,M,ped.phi2,Opt);
    case 'REMLHer'
        %Call NINGA_Her to run REML heritability analysis
        NINGA          = NINGA_Her(Y,M,ped.phi2,Opt);
    case 'Cov'
        %Call CovVec to run covariate inference
        NINGA          = CovVec(Y,M,C,ped.phi2,Opt);
    case 'GWAS'
       i          = str2num(getenv('SGE_TASK_ID'));
       i=1;
        %load pre-porcessed data, design matrix and kinship
       if ~Opt.PreProcess
           NINGA_PreProcess(Opt)
       end 
        load(fullfile(Opt.GWAo,sprintf('Misc_%s.mat',num2str(i))));       
        if ~isempty(Opt.Mask)
            Opt.Mask=Mask;
        end
        if Opt.inormal
            Y       = Y-repmat(mean(Y),size(Y,1),1);
            M       = M-repmat(mean(M),size(M,1),1);
            M (:,1) = ones(size(M,1),1);
            b       = M\Y;
            cov     = M*b;
            Y2       = inormal(Y-cov,'solar')+cov;
        end
        NINGA     = NINGA_GWA(Y ,M,i,Phi2,Opt);
end

%% calculate cluster FWE P-values and save the outputs
if ~strcmp(Opt.Method,'GWAS')
    %Write the outputs
    NINGA.Y.readwith  = Opt.Mask.readwith;
    NINGA.Y.extra     = Opt.Mask.extra;
    NINGA_write(NINGA,Opt);
    save(Opt.o,'NINGA','-v7.3')
else
    save(sprintf('%s_chr_%s.mat',Opt.o,num2str(i)),'NINGA','-v7.3')
    dlmwrite(sprintf('%s_chr_%s.txt',Opt.o,num2str(i)),NINGA.TS,' ')
end
end

