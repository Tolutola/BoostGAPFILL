function [RC] = FMconvert(train,libfmname)

datapath = strcat(pwd,'/data/');       %path of the data folder

%train = triu(train) - diag(diag(train));    %train is symmetric, delete redundant data
train = triu(train);

[r,c,v] = find(train);  %store row, col, value of train>0 data to libfmname, for further processing using python
RC = [r,c];
a = [v,r-1,c-1];
dlmwrite(strcat(datapath,libfmname),a,'delimiter',' ');
cd data;
cmd = sprintf('python processFM.py %s',libfmname);     %call python, convert to [v, r, c] to [v, r:1, c:1]
system(cmd);
cd ..;




