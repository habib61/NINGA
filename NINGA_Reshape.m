function [IndKin,IndData,idxData,idxGWA,idxKin,idxData2] = NINGA_Reshape(SubjID,KinID,GWAid)
%%This function recieves ID from kinship, data and gwa file, match then
%%returns indices. 
%%For future: modify this function to work only with uniqId and KinID 



%%
%Matching: remove subjects

KinID      = cellfun(@num2str,KinID,'UniformOutput',false);
GWAid      = cellfun(@num2str,GWAid,'UniformOutput',false);
SubjID     = cellfun(@num2str,SubjID,'UniformOutput',false);


delDataGWA = setxor(GWAid,SubjID); 

idxData    = false(size(SubjID));
idxGWA     = false(size(GWAid));
for s = 1:numel(delDataGWA),
    idxData = idxData | strcmpi(delDataGWA{s},SubjID);
    idxGWA  = idxGWA  | strcmpi(delDataGWA{s},GWAid);
end
SubjID(idxData)     = [];
GWAid(idxGWA)       = [];


delDataKin = setxor(GWAid,KinID);

idxKin     = false(size(KinID));
idxData2   = false(size(SubjID));

for s = 1:numel(delDataKin),
    idxKin  = idxKin  | strcmpi(delDataKin{s},KinID);
    idxData2= idxData2| strcmpi(delDataKin{s},SubjID);
end
KinID(idxKin)     = [];
SubjID(idxData2)  = [];
GWAid(idxData2)   = [];

%find rearranging indices with respect to gwaID
nS      = numel(GWAid);
IndKin  = zeros(nS,1);
IndData = zeros(nS,1);
for l=1:nS
    IndKin(l)   = find(cell2mat(cellfun(@(x) strcmpi(x,GWAid{l}), KinID,'UniformOutput',false)));
    IndData(l)  = find(cell2mat(cellfun(@(x) strcmpi(x,GWAid{l}), SubjID,'UniformOutput',false)));
end
end

