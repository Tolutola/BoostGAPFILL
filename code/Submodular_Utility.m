function benifit = Submodular_Utility(Ul,U,dA)

% for i=1:10
%     tmp = bsxfun(@minus,U,dA-Ul);
% end
tmp = bsxfun(@minus,U,dA-Ul);
tmp = min(tmp,0);
benifit = sum(tmp) - sum(min(Ul-dA,0));
