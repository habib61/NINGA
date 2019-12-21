function Opt = NINGA_defaults
% set up default parameters

%
%
% _____________________________________
% Habib Ganjgahi
% Statistic Department, uni of Oxford.
% Dec/2016


%% Set the input data

% text file contains address and name for inputdata
Opt.in            = 'dataS.txt';          
% Output string basename
%Opt.o             = 'LCope3CovEng';     
Opt.o             = 'Test4Eline';     

%desing and contrast file,
Opt.d             = 'Cov.csv';         % design matrix
Opt.c             = 'contEnq.csv';        % contrast file

%set the kinship matrix
Opt.kin           = 'phi2.csv';        % Kinship matrix
Opt.ped           = 'pedindex.csv';    % pedigree file for GRM it should be id file
Opt.subjID        = 'Subjid.csv';    % subject id 

              

% Analysis method, can be 'Cov' for fixed effect analysis or 'REMLHer' for
% heritability analysis
Opt.Method        = 'Cov';
 
% Number of permutations
Opt.nP            =  1000;                 

% Mask for image data, it could be empty
Opt.Mask.filename = [];            

%Calculate FWE corrected P-value for image elements
Opt.v.I           = true;    

%for spatial statistic set the below lines
%cluster forming threshold for cluster-wise analysis,
Opt.clus          = 2.3;%2.3;%chi2mixinv(.95,1);
%TFCE
Opt.tfce.I        = false;    %set this to "true" for TFCE
Opt.tfce.T2       = false;   %set this to "true" for TBSS or surface data

% surface and area file for performing spatial statistic of HCP data 
Opt.srf           = [];%'L.midthickness_MSMAll.32k_fs_LR.surf.gii';%'Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii';
Opt.srfare        = [];%'L_area.func.gii';

%for the GWA set the below lines:
% Adress for the folder that contains GWAS data
Opt.GWA           = ''; 
% Subjects Id in GWAS 
Opt.GWAid         = '';     

%% don't change below, these are advanced options
Opt.PreProcess    = false;
Opt.inormal       = false;

Opt.Permout       =  false;                % Write Permutation distribution
Opt.seed          =  100;                  % Seed for random number generator, in case of Meta analysis  
                                          

Opt.Perm.method   = 'Free';

if Opt.tfce.T2
    
    Opt.tfce.E       = 1;
    Opt.tfce.H       = 2;
    Opt.tfce.con     = 26;
    Opt.con          = 26;
else
    Opt.con          = 26;
    Opt.tfce.E       = 0.5;
    Opt.tfce.H       = 2;
    Opt.tfce.con     = 26;
end








