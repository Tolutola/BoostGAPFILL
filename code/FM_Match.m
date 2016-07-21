function [Lambda,scores] = FM_Match(hltrain, hltest,labels,k,ith_experiment,num_prediction)

[row,col] = size(hltrain);
[row,col2] = size(hltest);
datapath = strcat(pwd,'/data/');       %path of the data folder
FMtrainr = sprintf('FMtrainr_exp%d',ith_experiment);
%generate positive training
f1=fopen(strcat(datapath,FMtrainr),'w+');
for i = 1:col
    r = hltrain(:,i);
    ind = find(r);
    ind = [1,ind'-1];
    fprintf(f1,strcat(num2str(ind),'\r\n'));
end

for runthis = 1:0
%generate negative training (deleted, because performance is better if not generate negative)
load(strcat(datapath,'US1.mat'));   %load universal reactions which has been matched to S1's metabolites. 
for i = 1:col
    r = hltrain(:,i);
    num_meta = sum(r);
    while 1
        c = zeros(row,1);
        perm = randperm(size(hltrain,1));
        c(perm(1:num_meta)) = 1;
        %lia = ismember(c',US','rows');    %check if c in US
        lia = 0;      %to save time
        if lia == 0
            i
            break
        end
    end
    ind = [0,find(c)'-1];
    fprintf(f1,strcat(num2str(ind),'\r\n'));
end
fclose(f1);
end

FMtestr = sprintf('FMtestr_exp%d',ith_experiment);
%generate testing
f2=fopen(strcat(datapath,FMtestr),'w+');
for i = 1:col2
    r = hltest(:,i);
    ind = find(r);
    ind = [labels(i),ind'-1];
    fprintf(f2,strcat(num2str(ind),'\r\n'));
end
fclose(f2);
cd data;
cmd = sprintf('python processFM.py %s',FMtrainr);     
system(cmd);
cmd = sprintf('python processFM.py %s',FMtestr);
system(cmd);
cmd = sprintf('libfm.exe -task c -train %s.libfm -test %s.libfm -dim "1,1,%d" -out outr_%d.txt -iter 1000',FMtrainr,FMtestr,k,ith_experiment);
evalc('system(cmd);');    %run libFM in silence mode
pred = load(sprintf('outr_%d.txt',ith_experiment));    %load the output file of libFM
delete(FMtrainr);
delete(FMtestr);
delete(sprintf('outr_%d.txt',ith_experiment));
delete(sprintf('%s.libfm',FMtrainr));
delete(sprintf('%s.libfm',FMtestr));
cd ..;

scores = pred;
Lambda = zeros(col2,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end
