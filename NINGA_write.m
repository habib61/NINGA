function NINGA_write(NINGA,Opt)
%Should recieve NINGA and Opt and writes corrected and uncorrected Pvalues
%and parmater estimates and their variances, 

%% Calculate FWE corrected P-values if requested
% for c=1:Opt.nCon 
%     %Cluster-wise FWE P-value
%     if ~isempty(Opt.clus)
%         if isempty(NINGA.Clus.Size{c})
%             NINGA.Clus.FWEp     = 0;
%             NINGA.Clus.ImgP{c}  = zeros(size(Opt.Mask.data));
%         else
%             for j=1:numel(NINGA.Clus.Size{c})
%                 NINGA.Clus.FWEp{c}(1,j) = mean(NINGA.Clus.Size{c}(j)>=NINGA.Clus.Max{c});
%             end
%             NINGA.Clus.ImgP{c}                  = zeros(size(Opt.Mask.data));
%             NINGA.Clus.ImgP{c}(Opt.Mask.data>0) = NINGA.Clus.FWEp{c};
%         end      
%     end        
% end
     
%% 
readwith = NINGA.Y.readwith;
extra    = NINGA.Y.extra;
if strcmp(readwith,'load')
   Opt.Mask.data=[];
end


for c=1:Opt.nCon
    if strcmp(Opt.Method,'Cov')   
        S          = NINGA_maskstruct(NINGA.Beta{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_',num2str(c),'_','Beta');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.BetaVar{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_',num2str(c),'_','BetaVar');
        NINGA_miscwrite(S)
        
        S          = NINGA_maskstruct(NINGA.TS{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_',num2str(c),'_','TS');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.P.Para{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_',num2str(c),'_','Uncor_Para');
        NINGA_miscwrite(S)

        if Opt.nP>0
            S          = NINGA_maskstruct(NINGA.P.Uncor{c},Opt.Mask.data,readwith,extra);
            S.filename = horzcat(Opt.o,'_',num2str(c),'_','Uncor_Perm');
            NINGA_miscwrite(S)
            
            if Opt.v.I
                S          = NINGA_maskstruct(NINGA.P.FWE{c},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_',num2str(c),'_','FWE_P');
                NINGA_miscwrite(S)
            end
            %Cluster
            if ~isempty(Opt.clus)
                S          = NINGA_maskstruct(NINGA.Clus.FWEp{c},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_',num2str(c),'_','c');
                NINGA_miscwrite(S)
            end
            %TFCE
            if Opt.tfce.I
                S          = NINGA_maskstruct(NINGA.TFCE.FWE{1},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','tfce_fwep');
                NINGA_miscwrite(S)
                
                S          = NINGA_maskstruct(NINGA.TFCE.FWE{1},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','tfce_uncorp');
                NINGA_miscwrite(S)
            end 
        end        
    else            
        S          = NINGA_maskstruct(NINGA.h2{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_','h2');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.h2Var{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_','h2Var');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.TS{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_','TS');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.SigmaY{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_','SigmaP');
        NINGA_miscwrite(S)

        S          = NINGA_maskstruct(NINGA.P.Para{c},Opt.Mask.data,readwith,extra);
        S.filename = horzcat(Opt.o,'_','Uncor_Para');
        NINGA_miscwrite(S)

        if Opt.nP>0
            S          = NINGA_maskstruct(NINGA.P.Uncor{c},Opt.Mask.data,readwith,extra);
            S.filename = horzcat(Opt.o,'_','Uncor_Perm');
            NINGA_miscwrite(S)
            %Voxel or element-wise FWE
            if Opt.v.I
                S          = NINGA_maskstruct(NINGA.P.FWE{c},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','FWE_P');
                NINGA_miscwrite(S)
            end
            %Cluster
            if ~isempty(Opt.clus)
                S          = NINGA_maskstruct(NINGA.Clus.FWEp{1},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','c');
                NINGA_miscwrite(S)
            end
            %TFCE
            if Opt.tfce.I
                S          = NINGA_maskstruct(NINGA.TFCE.FWE{1},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','tfce_fwep');
                NINGA_miscwrite(S)
                
                S          = NINGA_maskstruct(NINGA.TFCE.FWE{1},Opt.Mask.data,readwith,extra);
                S.filename = horzcat(Opt.o,'_','tfce_uncorp');
                NINGA_miscwrite(S)
            end
        end
    end        
end

% this bit is for writing cov sigmaY, check it out when testing cov
% analysis
%             if ~isempty(Opt.Mask.data)
%                 dim  = [size(Opt.Mask.data) size(NINGA.SigmaY{c},1)];
%                 kk   = zeros(dim(1:3));
%                 vol  = zeros(dim);
%                 for i=1:dim(4)
%                     kk(Opt.Mask.data>0) = NINGA.SigmaY{c}(i,:);
%                     vol(:,:,:,i)         = kk;
%                 end
%                 S         = NINGA_maskstruct(vol,Opt.Mask.data,readwith,extra);
% 
%             else
%                 S          = NINGA_maskstruct(NINGA.SigmaY{c},Opt.Mask.data,readwith,extra);
%             end
end