function [Lambda,scores] = HPLSF(wtrain,test,k,ith_experiment,hltrain,hltest,num_prediction)


%evalc('FM(wtrain,test,k,ith_experiment,1);');

[Feature,trainlabels,Feature1] = readFMmodel(hltrain,hltest,wtrain,k);


X = [Feature;Feature1];
Z = zscore(X);
Z(:,1:size(Z,2)-rank(Z)) = [];
Feature = Z(1:size(Feature,1),:)+eps;
Feature1 = Z(size(Feature,1)+1:end,:)+eps;

lr = mnrfit(Feature,trainlabels+1);
scores = mnrval(lr,Feature1);
scores(:,1) = [];

Lambda = zeros(size(hltest,2),1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end


