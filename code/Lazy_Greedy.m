function [benifit,order,nextselect] = Lazy_Greedy(Ul,U,dA,order,util)

util1 = min(Ul+U(:,order(1))-dA,0);
util1 = sum(util1) - sum(min(Ul-dA,0));

if util1>=util(2)
    nextselect = order(1);
    order = order(2:end);
    benifit = util(2:end);
else
    tmp = Submodular_Utility(Ul,U,dA);
    [benifit,order] = sort(tmp,2,'descend');
    nextselect = order(1);
    order = order(2:end);
    benifit = benifit(2:end);
end