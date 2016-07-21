function [Lambda,scores] = HCN(wtrain,test,hltest,num_prediction)

[~,cc] = size(hltest);
train = spones(wtrain);
sim = train*train;
sim = sim - diag(diag(sim));
sim = triu(sim);
scores = zeros(cc,1);
for i = 1:cc
    r = hltest(:,i);
    rA = r*r';
    rA = rA - diag(diag(rA));
    scores(i) = sum(sum(rA.*sim))/(nnz(r)*(nnz(r)-1)/2);
end

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end