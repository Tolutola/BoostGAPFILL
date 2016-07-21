%  This program is for preprocessing the universal reactions for a specific metabolic model. It also unifies the data formats between universal reactions and the specific model.
%  You should first run the Meta_Download_U.m
clear all;

%% load the current model
datapath = strcat(pwd,'/data/');

% select a model to read, you can add your models here. For each new model, only run this one time.
mo = 9;

for mo = 1:9

switch mo
    case 1
        model = 'iJO1366';
    case 2
        model = 'iAF1260b';
    case 3
        model = 'e_coli_core';
    case 4
        model = 'RECON1';
    case 5
        model = 'iAT_PLT_636';
    case 6
        model = 'iAB_RBC_283';
    case 7
        model = 'iAF692';
    case 8
        model = 'iHN637';
    case 9
        model = 'iIT341';
end

load(sprintf('data/%s.mat',model),model);  [~,Model] = evalc(model);

S = Model.S;    %stoichiometric matrix S
S = full(spones(S));   %make to binary
R = Model.rxns;   %reactions R
M = Model.mets;   %metabolites M

%% generate universal S matrix US for the current Model

load(strcat(datapath,'unrmet.mat'));   %load the universal reactions' metabolites
US = zeros(size(M,1),length(unrmet));    %universal reactions' S matrix (ONLY KEEP METABOLITES that appear IN ORIGINAL M)
num_met_ur = [];     %number of metabolites in each unr
for i = 1:length(unrmet)
    i;
    cr = unrmet{i};     %current reaction
    cmet = cr.metabolites;           %current metabolites
    dS = zeros(size(M));
    for j = 1:length(cmet)    %for each metabolite of the reaction, compare with M and keep the matches and form dS
        cname = strcat(cmet{j}.bigg_id,'_',cmet{j}.compartment_bigg_id);     %unify the format
        dS = dS + strcmp(M,cname)*(cmet{j}.stoichiometry);
    end
    num_met_ur = [num_met_ur,length(cmet)];    %record the number of metabolites for each reaction in unr
    US(:,i) = dS;    %append dS to US
end
LUS = logical(US);
ind_US = find(sum(LUS)==num_met_ur);    %only keep those reactions (columns) in US whose metabolites are all in M
US = US(:,ind_US);         %the reduced US
tmp = importdata(strcat(datapath,'unrnames.txt'))';    %import universal reactions' names
unrnames = tmp(1,ind_US);   %the reduced unr names
Model.US = US;
Model.unrnames = unrnames;

%% for regularize Model reactions' names, make them the same format with universal reactions
for i=1:size(R,1)
    rname = R{i};
    if length(rname)<5
        continue
    end
    if rname(end-4:end-1)=='copy'      %delete the "copy" from the reaction id
        rname = rname(1:end-6);
        R{i} = rname;
    end
end
Model.rxns = R;

%% check if all R have corresponding names in unr
tmp = importdata(strcat(datapath,'unrnames.txt'))';    %import universal reactions' names
a = strcmp(repmat(R,1,size(tmp,2)),repmat(tmp,size(R,1),1));
b = sum(a,2);
nomatch = R(b==0);
assert(isempty(nomatch));

%% save all needed file into the original model file
save(sprintf('data/%s.mat',model),'Model');

end
