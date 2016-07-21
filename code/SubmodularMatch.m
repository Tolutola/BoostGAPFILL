function [Lambda,scores] = SubmodularMatch(wtrain,test,k,ith_experiment,hltest,num_prediction)

[~,~,~,dA] = evalc('FM(wtrain,test,k,ith_experiment);');
mask = +(spones(wtrain)==0);
mask = mask(:);
mask = sparse(mask);
 
% mask = ~logical(spones(wtrain));
% mask = mask(:);

[rr,cc] = size(hltest);
U = [];
for i=1:cc
    u = hltest(:,i);
    u = u*u';
    u = sparse(u(:));    
    U = [U,u];
end

dA = dA(:);
[Y,~] = sort(dA,1,'descend');
select = round(0.1*size(dA,1));
dA(dA<Y(select)) = 0;
dA = bsxfun(@times,dA,mask);
U = bsxfun(@times,U,mask);
dA1 = dA;
Lambda = zeros(cc,1);
Ul = sparse(zeros(rr^2,1));
for i = 1:num_prediction
    i
    if i==1
        benifit = Submodular_Utility(Ul,U,dA1);
        [benifit,order] = sort(benifit,2,'descend');
        ind = order(1);
        order = order(2:end);
        benifit = benifit(2:end);
        Lambda(ind) = 1;
    else  
        %U(:,ind) = -1000*ones(rr^2,1);
        U(1,ind) = -500;
        [benifit,order,ind] = Lazy_Greedy(Ul,U,dA1,order,benifit);
        Lambda(ind) = 1;
    end
    Ul = Ul + U(:,ind);
end
scores = Lambda;
Lambda = logical(Lambda);

