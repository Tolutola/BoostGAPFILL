function [Lambda,scores] = ILSQ_Alternate(wtrain,test,k,hltest,num_prediction,max_iter,blacklist,solver)
%  wtrain: the weighted train adjacency matrix
%  test: no use, for consistency with other link prediction algorithms (can be ignored)
%  k: the number of latent factors used in matrix factorization
%  max_iter: the maximum number of iterations allowed
%%
mask = spones(wtrain)==0;
[rr,cc] = size(hltest);


U1 = [];
for i=1:cc
    u = hltest(:,i);
    u = u*u';
    u = sparse(u(:));
    U1 = [U1,u];
end


scores = zeros(cc,1);   %initialize scores
DS = zeros(1,max_iter);  %record the ||newscores - scores||
Scores = [];

for iter = 1:max_iter
    
    
    %completion step
    wtrain1 = wtrain + (hltest*diag(scores)*hltest').*(~mask);
    ith_experiment=1; % used to prevent file conflict
    [~,~,~,wdA] = evalc('FM(wtrain1,test,k,ith_experiment);');
    
    %matching step
    [~,newscores] = ILSQ_Match(wdA,U1,num_prediction,rr,cc,mask,0,blacklist,solver);
    
    %stop criterion
    ds = norm(newscores - scores);
    DS(iter) = ds;
    if ds>min(DS(1:max(iter-1,1)))        %if stops because of stop criterion, calculate the mean scores before this iteration
        break
    end
    scores = newscores;
    Scores = [Scores,scores];
end
scores = mean(Scores(:,1:max(end-1,1)),2);

Lambda = zeros(cc,1);
[~,I] = sort(scores,1,'descend');
if num_prediction > 1
    if num_prediction> length(I), num_prediction=length(I); end
    Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
    Lambda = logical(Lambda);
else
    Lambda(scores>num_prediction) = 1;
    Lambda = logical(Lambda);
end

