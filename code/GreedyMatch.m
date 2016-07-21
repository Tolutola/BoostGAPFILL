function [Lambda,scores] = GreedyMatch(wtrain,test,k,ith_experiment,hltest,num_prediction)

[~,~,dA,~] = evalc('FM(wtrain,test,k,ith_experiment);');

[rr,cc] = size(hltest);
scores = zeros(cc,1);
for i=1:cc
    i;
    r = hltest(:,i);
    rA = r*r';
    %rA = rA - diag(diag(rA));
    %Lambda(i) = (nnz(rA)-nnz(rA.*dA))/nnz(rA)<=0.3;
    scores(i) = (nnz(rA.*dA))/nnz(rA);
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


