
function ped = NINGA_kin(pedfile,phifile)
% Load the pedigree
% 
% Usage for pedigree:
% pedigree = loadped('pedindex.csv','phi2.csv')
% 
% Usage for empiricalkinship:
% pedigree = loadped('phi2.csv')

% The input files must have been converted to CSV files.
% The output is a struct that contains the headers for
% the pedindex, the pedindex itself (as a double array)
% and the matrices Phi2 and Delta7.
% 
% _____________________________________
% Anderson M. Winkler., Habib Ganjgahi
% The University of Warwick/ Statistic Department
% Jan 2015

[pathstr,name,ext] = fileparts(pedfile);
if strcmpi(ext,'.csv') %it means that the phi2 and pedigree are in Solar format
    % Read the pedindex.csv file
    fid = fopen(pedfile,'r');
    pedindex = textscan(fid,repmat('%s',[1 8]),'Delimiter',',');
    fclose(fid);
    % Move the headers to a separate variable
    headers = cell(size(pedindex'));
    for c = 1:numel(pedindex),
        headers{c}     = pedindex{c}{1};
        pedindex{c}(1) = [];
    end
    % Subject IDs
    cid = strcmpi('ID',headers);
    id  = pedindex{cid};

    % Convert char to num and to an array
    pidx = zeros(numel(pedindex{1}),numel(pedindex));
    for c = 1:numel(pedindex),
        pedindex{c} = str2double(pedindex{c}(:));
        pidx(:,c) = pedindex{c}(:);
    end

    % Read the phi2.csv file
    fid      = fopen(phifile,'r');
    phitable = textscan(fid,repmat('%f',[1 4]),'Delimiter',',');
    fclose(fid);

    % Organise as simple table
    phitable = [phitable{1} phitable{2} phitable{3} phitable{4}];

    % Reorganise as symmetric matrices
    tmpvar = phitable(:,1:2); N = max(tmpvar(:));
    Phi2 = zeros(N);    % Kinship
    Delta7 = zeros(N);  % Jacquard D7 coefficients

    % For each row in the phi2.gz file
    for r = 2:size(phitable,1),
        Phi2(phitable(r,1),phitable(r,2)) = phitable(r,3);
        Phi2(phitable(r,2),phitable(r,1)) = phitable(r,3);
        Delta7(phitable(r,1),phitable(r,2)) = phitable(r,4);
        Delta7(phitable(r,2),phitable(r,1)) = phitable(r,4);
    end
    % Organize the output
    ped.hdr  = headers;
    ped.id   = id;
    ped.pidx = pidx;
    ped.phi2 = Phi2;
    ped.d7   = Delta7;
else
%     %grm from binary file
%     tmp = dlmread(pedfile);
%     pedindex = tmp(:,1);
%     Phi2 = dlmread(phifile);
%     ped.phi2 = Phi2;
%     ped.id   = pedindex;
    fid      = fopen(pedfile,'r');
    pedindex = textscan(fid,'%s %s');
    fclose(fid);
    
    % Read the phi2.csv file
%    phifile  = gunzip(fullfile(pwd,phifile));
    fid      = fopen(phifile,'r');
    phitable = textscan(fid,'%f %f %s %f');
    fclose(fid);
    
    phitable          = [phitable{1} phitable{2} phitable{4}];

    tmpvar   = phitable(:,1:2); N = max(tmpvar(:));
    Phi2     = zeros(N);    % Kinship
    % For each row in the phi2.gz file
    for r = 1:size(phitable,1),
        Phi2(phitable(r,1),phitable(r,2)) = phitable(r,3);
        Phi2(phitable(r,2),phitable(r,1)) = phitable(r,3);
    end

    ped.phi2 = Phi2;
    ped.id   = pedindex{2};
end
end
