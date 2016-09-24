function new_metabolites= newMetFind(newReactions,M,unrnames)
%% MUnames,
load('data/unrmet.mat');   %load the universal reactions' metabolites
new_metabolites = {};   %to record new metabolites
exotic_number = [];   %to record the number of exotic metabolites in each selected reaction
cellind = 1;
M1 = M';
for sr = 1:length(newReactions)
    srname = newReactions{sr};
    usr = strcmp(srname,unrnames);
    cr = unrmet{usr};
    cmet = cr.metabolites;           %current metabolites
    num_exist = 0;
    for j = 1:length(cmet)    %for each metabolite of the reaction, compare with M and keep the matches and form dS
        cname = strcat(cmet{j}.bigg_id,'_',cmet{j}.compartment_bigg_id);     %unify the format
        tmp = strcmp(cname,M1);
        tmp = sum(tmp);
        num_exist = num_exist + tmp;
        if tmp==0
            new_metabolites{cellind} = cname;
            cellind = cellind + 1;
        end
    end
    M1 = [M',new_metabolites];
    exotic_number = [exotic_number,length(cmet)-num_exist];
    
end
